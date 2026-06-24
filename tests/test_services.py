"""
Integration tests for external service endpoints.

Run with:  pytest tests/test_services.py -v
These make real HTTP requests to detect outages and API changes.
"""

import pytest
import requests

# ---------------------------------------------------------------------------
# Stable dataset IDs used as fixtures across tests
# ---------------------------------------------------------------------------
MASSIVE_ID = "MSV000081351"
METABOLIGHTS_ID = "MTBLS1"
WORKBENCH_ID = "ST000001"
FASST_USI = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00006582001"
FASST_HOST = "https://api.fasst.gnps2.org"
TIMEOUT = 20


# ---------------------------------------------------------------------------
# MassIVE (PROXI)
# ---------------------------------------------------------------------------

class TestMassiveService:
    URL = f"https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/{MASSIVE_ID}"

    def test_endpoint_reachable(self):
        r = requests.get(self.URL, timeout=TIMEOUT)
        assert r.status_code == 200, f"Expected 200, got {r.status_code}"

    def test_response_is_json(self):
        r = requests.get(self.URL, timeout=TIMEOUT)
        data = r.json()
        assert isinstance(data, dict)

    def test_response_has_expected_keys(self):
        r = requests.get(self.URL, timeout=TIMEOUT)
        data = r.json()
        for key in ("title", "summary", "species", "publications"):
            assert key in data, f"Missing key: '{key}'"

    def test_title_is_nonempty_string(self):
        r = requests.get(self.URL, timeout=TIMEOUT)
        title = r.json().get("title", "")
        assert isinstance(title, str) and len(title) > 0


# ---------------------------------------------------------------------------
# MetaboLights
# ---------------------------------------------------------------------------

class TestMetaboLightsService:
    URL = f"https://www.ebi.ac.uk:443/metabolights/ws/studies/{METABOLIGHTS_ID}?investigation_only=true"

    def test_endpoint_reachable(self):
        r = requests.get(self.URL, timeout=TIMEOUT)
        assert r.status_code == 200, f"Expected 200, got {r.status_code}"

    def test_response_is_json(self):
        r = requests.get(self.URL, timeout=TIMEOUT)
        data = r.json()
        assert isinstance(data, dict)

    def test_response_has_expected_keys(self):
        r = requests.get(self.URL, timeout=TIMEOUT)
        data = r.json()
        assert "isaInvestigation" in data, "Missing key: 'isaInvestigation'"
        assert "mtblsStudy" in data, "Missing key: 'mtblsStudy'"

    def test_study_list_nonempty(self):
        r = requests.get(self.URL, timeout=TIMEOUT)
        studies = r.json().get("isaInvestigation", {}).get("studies", [])
        assert len(studies) > 0, "Expected at least one study in response"


# ---------------------------------------------------------------------------
# Metabolomics Workbench
# ---------------------------------------------------------------------------

class TestWorkbenchService:
    URL = f"https://www.metabolomicsworkbench.org/rest/study/study_id/{WORKBENCH_ID}/summary"

    def test_endpoint_reachable(self):
        r = requests.get(self.URL, timeout=TIMEOUT)
        assert r.status_code == 200, f"Expected 200, got {r.status_code}"

    def test_response_is_json(self):
        r = requests.get(self.URL, timeout=TIMEOUT)
        data = r.json()
        assert isinstance(data, dict)

    def test_response_has_expected_keys(self):
        r = requests.get(self.URL, timeout=TIMEOUT)
        data = r.json()
        for key in ("study_id", "study_title", "species"):
            assert key in data, f"Missing key: '{key}'"

    def test_study_id_matches(self):
        r = requests.get(self.URL, timeout=TIMEOUT)
        assert r.json().get("study_id") == WORKBENCH_ID


# ---------------------------------------------------------------------------
#  ReDU metadata dump / FASST / GNPS2
# ---------------------------------------------------------------------------

REDU_EXPECTED_COLUMNS = {
    "filename",
    "ATTRIBUTE_DatasetAccession",
    "NCBITaxonomy",
    "UBERONBodyPartName",
    "DOIDCommonName",
    "BiologicalSex",
    "AgeInYears",
}

class TestReDUService:
    URL = "https://redu.gnps2.org/dump"

    def test_endpoint_reachable(self):
        with requests.get(self.URL, stream=True, timeout=TIMEOUT) as r:
            assert r.status_code == 200, f"Expected 200, got {r.status_code}"

    def test_response_is_nonempty(self):
        with requests.get(self.URL, stream=True, timeout=TIMEOUT) as r:
            first_chunk = next(r.iter_content(chunk_size=1024), b"")
        assert len(first_chunk) > 0, "Response body is empty"

    def test_header_row_has_expected_columns(self):
        """Read only the first line to check TSV columns without downloading the full file."""
        with requests.get(self.URL, stream=True, timeout=TIMEOUT) as r:
            r.raise_for_status()
            first_line = b""
            for chunk in r.iter_content(chunk_size=4096):
                first_line += chunk
                if b"\n" in first_line:
                    break

        header_line = first_line.split(b"\n")[0].decode("utf-8")
        columns = set(header_line.strip().split("\t"))
        missing = REDU_EXPECTED_COLUMNS - columns
        assert not missing, f"Missing expected columns: {missing}"


# ---------------------------------------------------------------------------
class TestFASSTService:

    def test_search_endpoint_accepts_job(self):
        """POST to /search should return a task ID immediately."""
        params = {
            "library": "gnpsdata_index",
            "usi": FASST_USI,
            "analog": "No",
            "pm_tolerance": 0.02,
            "fragment_tolerance": 0.02,
            "cosine_threshold": 0.7,
        }
        r = requests.post(f"{FASST_HOST}/search", json=params, timeout=TIMEOUT)
        assert r.status_code == 200, f"Expected 200, got {r.status_code}"
        data = r.json()
        assert "id" in data, f"Response missing 'id' key: {data}"

    def test_result_endpoint_reachable(self):
        """Submit a job and confirm the result endpoint returns a valid status."""
        params = {
            "library": "gnpsdata_index",
            "usi": FASST_USI,
            "analog": "No",
            "pm_tolerance": 0.02,
            "fragment_tolerance": 0.02,
            "cosine_threshold": 0.7,
        }
        r = requests.post(f"{FASST_HOST}/search", json=params, timeout=TIMEOUT)
        task_id = r.json()["id"]

        r2 = requests.get(f"{FASST_HOST}/search/result/{task_id}", timeout=TIMEOUT)
        assert r2.status_code == 200, f"Result endpoint returned {r2.status_code}"
        data = r2.json()
        # Either still pending or already done — both are valid shapes
        assert isinstance(data, dict), "Expected dict response from result endpoint"
        if data.get("status") == "PENDING":
            pass  # service is up and processing
        else:
            assert "results" in data, f"Completed response missing 'results' key: {data}"
