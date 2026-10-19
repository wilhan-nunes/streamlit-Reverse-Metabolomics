"""
Per-service reporting for tests/test_services.py.

Each test class declares the external service it checks (SERVICE) and whether
the app depends on it critically (CRITICAL). Failures are grouped by service
and summarised at the end of the run, so a broken endpoint shows up as
"MassIVE: SSL certificate verification failed" instead of a requests traceback.
On GitHub Actions the summary is also written to the job page and each failing
service gets its own error annotation.
"""

import os
import re

import pytest
import requests

# service -> {"critical": bool, "passed": int, "failed": [test names], "reasons": [str]}
_results = {}


def _describe(exc):
    """Reduce an exception from a service check to a one-line reason."""
    errors = requests.exceptions
    if isinstance(exc, errors.SSLError):
        match = re.search(r"verify failed: (.+?) \(_ssl", str(exc))
        detail = f": {match.group(1)}" if match else ""
        return f"SSL certificate verification failed{detail}"
    if isinstance(exc, errors.Timeout):
        return "request timed out"
    if isinstance(exc, errors.ConnectionError):
        return "could not connect (DNS failure, connection refused or reset)"
    if isinstance(exc, errors.JSONDecodeError):
        return "response was not valid JSON"
    if isinstance(exc, AssertionError):
        return str(exc).split("\n")[0] or "assertion failed"
    return f"{type(exc).__name__}: {exc}"


@pytest.hookimpl(wrapper=True)
def pytest_runtest_makereport(item, call):
    report = yield
    service = getattr(item.cls, "SERVICE", None)
    if service is None or report.when != "call":
        return report

    entry = _results.setdefault(service, {
        "critical": getattr(item.cls, "CRITICAL", False),
        "passed": 0,
        "failed": [],
        "reasons": [],
    })
    if report.passed:
        entry["passed"] += 1
    elif report.failed:
        entry["failed"].append(item.name)
        reason = _describe(call.excinfo.value)
        if reason not in entry["reasons"]:
            entry["reasons"].append(reason)
    return report


def _annotation_escape(text):
    return text.replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")


def pytest_terminal_summary(terminalreporter):
    if not _results:
        return
    tr = terminalreporter
    # Critical services first, then alphabetical
    ordered = sorted(_results.items(), key=lambda kv: (not kv[1]["critical"], kv[0]))
    failing = [(name, e) for name, e in ordered if e["failed"]]
    critical_failing = [name for name, e in failing if e["critical"]]

    tr.section("external service status")
    if critical_failing:
        tr.write_line(f"CRITICAL services failing: {', '.join(critical_failing)}", red=True, bold=True)
    elif failing:
        tr.write_line("Critical services (FASST, ReDU) OK; non-critical services failing.", yellow=True, bold=True)

    md = ["| Service | Status | Reason |", "|---|---|---|"]
    for name, e in ordered:
        label = f"{name} [CRITICAL]" if e["critical"] else name
        total = e["passed"] + len(e["failed"])
        if e["failed"]:
            reason = "; ".join(e["reasons"])
            tr.write_line(f"FAIL  {label}: {len(e['failed'])}/{total} checks failed - {reason}", red=True)
            md.append(f"| {label} | ❌ {len(e['failed'])}/{total} failed | {reason} |")
        else:
            tr.write_line(f"OK    {label}: {total}/{total} checks passed", green=True)
            md.append(f"| {label} | ✅ OK | |")

    if os.environ.get("GITHUB_ACTIONS") == "true":
        for name, e in failing:
            title = f"{name} is failing" + (" (CRITICAL)" if e["critical"] else "")
            message = f"{len(e['failed'])} check(s) failed: " + "; ".join(e["reasons"])
            tr.write_line(f"::error title={_annotation_escape(title)}::{_annotation_escape(message)}")

    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary_path:
        heading = "## External service status\n\n"
        if critical_failing:
            heading += f"**🚨 Critical services failing: {', '.join(critical_failing)}**\n\n"
        with open(summary_path, "a", encoding="utf-8") as fh:
            fh.write(heading + "\n".join(md) + "\n")
