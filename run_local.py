#!/usr/bin/env python3
"""Launch the Reverse Metabolomics app on this machine.

Usage:  uv run run_local.py       (or: python run_local.py)

If the dependencies are not installed, this script re-launches itself under
`uv`, which installs a private copy of Python and every requirement into a
cache and leaves the rest of the system untouched.

This is the local counterpart to run_server.sh. The differences that matter:
the server binds 0.0.0.0 on a fixed port behind nginx, this binds 127.0.0.1 on
a free port and opens your browser at it.
"""

import importlib.util
import os
import shutil
import socket
import subprocess
import sys
import threading
import time
import webbrowser
from pathlib import Path

REPO_DIR = Path(__file__).resolve().parent
BOOTSTRAP_FLAG = "REVMET_BOOTSTRAPPED"
PYTHON_VERSION = "3.12"


def data_dir() -> Path:
    """Where the ReDU metadata file and other downloads are kept.

    Next to the user's home rather than the working directory, so the ~1 GB
    download survives moving or re-cloning the repo.
    """
    override = os.environ.get("REVMET_DATA_DIR")
    directory = Path(override) if override else Path.home() / ".reverse-metabolomics"
    directory.mkdir(parents=True, exist_ok=True)
    return directory


# Every third-party module the app imports at startup. Checking only one of
# these is not enough: an unrelated virtualenv that happens to have Streamlit
# would satisfy it, and the app would then fail on the first missing import.
REQUIRED_MODULES = (
    "streamlit", "pandas", "numpy", "scipy", "seaborn",
    "matplotlib", "requests", "tqdm", "pyarrow",
)


def missing_modules() -> list:
    """Names of required modules that cannot be found in this environment."""
    absent = []
    for name in REQUIRED_MODULES:
        try:
            if importlib.util.find_spec(name) is None:
                absent.append(name)
        except (ImportError, ValueError):
            absent.append(name)
    return absent


def bootstrap_with_uv(absent) -> int:
    """Re-exec under `uv run`, which resolves Python and the requirements."""
    uv = shutil.which("uv")
    if uv is None:
        sys.exit(
            "Dependencies are not installed, and `uv` was not found.\n\n"
            "Install uv, then run this script again:\n"
            "  macOS/Linux:  curl -LsSf https://astral.sh/uv/install.sh | sh\n"
            "  Windows:      powershell -c \"irm https://astral.sh/uv/install.ps1 | iex\"\n\n"
            "Or install the requirements into your own environment:\n"
            f"  pip install -r {REPO_DIR / 'requirements.txt'}"
        )

    # flush: this is followed by a slow subprocess, and a buffered notice would
    # reach the user only after the wait it is meant to explain.
    print(f"Missing dependencies: {', '.join(absent)}", flush=True)
    print("Setting up (first run downloads Python and the dependencies)...", flush=True)

    # Build the environment from requirements.txt alone. An active virtualenv or
    # conda env would otherwise be reused, which is how we got here.
    env = dict(os.environ, **{BOOTSTRAP_FLAG: "1"})
    env.pop("VIRTUAL_ENV", None)
    env.pop("CONDA_PREFIX", None)
    env.pop("PYTHONPATH", None)

    return subprocess.call(
        [
            uv, "run",
            "--no-project",
            "--python", PYTHON_VERSION,
            "--with-requirements", str(REPO_DIR / "requirements.txt"),
            "python", str(REPO_DIR / "run_local.py"),
        ],
        env=env,
    )


def free_port() -> int:
    """Ask the OS for an unused port so two copies can run side by side."""
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def open_browser_when_ready(port: int, timeout: float = 60.0) -> None:
    """Wait for the server to accept connections, then open the browser.

    Streamlit only opens the browser itself in non-headless mode, which also
    triggers its interactive "enter your email" prompt on first run. We stay
    headless to skip that and do this part ourselves.
    """
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=0.5):
                break
        except OSError:
            time.sleep(0.25)
    else:
        return  # Never came up; the URL is already printed as a fallback.

    try:
        webbrowser.open(f"http://127.0.0.1:{port}")
    except Exception:
        pass  # Non-fatal: the user can click the printed URL.


def main() -> int:
    absent = missing_modules()
    if absent:
        if os.environ.get(BOOTSTRAP_FLAG):
            sys.exit(
                f"Still missing after setup: {', '.join(absent)}\n"
                "Please report this, including the output above, at\n"
                "https://github.com/wilhan-nunes/streamlit-Reverse-Metabolomics/issues"
            )
        return bootstrap_with_uv(absent)

    # Run from the repo regardless of where the user invoked us, so relative
    # lookups (the git hash in the About menu, example_data/) resolve the same
    # way whether this was launched from a terminal or by double-clicking.
    os.chdir(REPO_DIR)

    os.environ.setdefault("REVMET_DATA_DIR", str(data_dir()))
    os.environ["REVMET_LOCAL"] = "1"

    port = free_port()
    print(f"\nStarting Reverse Metabolomics at http://127.0.0.1:{port}")
    print("Your browser should open automatically.")
    print("Press Ctrl+C in this window to stop.\n", flush=True)

    threading.Thread(target=open_browser_when_ready, args=(port,), daemon=True).start()

    from streamlit.web import cli as stcli

    sys.argv = [
        "streamlit", "run", str(REPO_DIR / "app.py"),
        "--server.address", "127.0.0.1",
        "--server.port", str(port),
        "--server.headless", "true",
        "--browser.gatherUsageStats", "false",
    ]
    try:
        return stcli.main()
    except KeyboardInterrupt:
        print("\nStopped.")
        return 0


if __name__ == "__main__":
    sys.exit(main())
