import streamlit as st
import requests
from typing import Dict, Any, List, Optional
from enum import Enum


class DataRepository(Enum):
    """Enum for supported data repositories"""

    MASSIVE = "massive"
    METABOLIGHTS = "metabolights"
    WORKBENCH = "workbench"


def detect_repository(dataset_id: str) -> Optional[DataRepository]:
    """
    Detect which repository a dataset ID belongs to.

    Parameters:
    -----------
    dataset_id : str
        The dataset identifier

    Returns:
    --------
    DataRepository or None
        The detected repository or None if not recognized
    """
    dataset_id_upper = dataset_id.upper()

    if dataset_id_upper.startswith("MSV"):
        return DataRepository.MASSIVE
    elif dataset_id_upper.startswith("MTBLS"):
        return DataRepository.METABOLIGHTS
    elif dataset_id_upper.startswith("ST"):
        return DataRepository.WORKBENCH
    else:
        return None


def fetch_massive_metadata(msv_id: str) -> Dict[str, Any]:
    """Fetch metadata from MassIVE repository"""
    url = f"http://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/{msv_id}"
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    return response.json()


def fetch_metabolights_metadata(mtbls_id: str) -> Dict[str, Any]:
    """Fetch metadata from MetaboLights repository"""
    url = f"https://www.ebi.ac.uk:443/metabolights/ws/studies/{mtbls_id}?investigation_only=true"
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    return response.json()


def fetch_workbench_metadata(st_id: str, summary_only: bool = False) -> Dict[str, Any]:
    """Fetch metadata from Metabolomics Workbench repository"""
    url = f"https://www.metabolomicsworkbench.org/rest/study/study_id/{st_id}/"
    if summary_only:
        url += "summary"
    else:
        url += "mwtab"
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    return response.json()


def display_massive_metadata(data: Dict[str, Any], dataset_id: str) -> None:
    """Display MassIVE metadata"""

    # Title
    st.header(data.get("title", "No title available"))

    # Dataset ID badge
    st.markdown(f"**Dataset ID:** `{dataset_id}`")

    # Summary
    st.subheader("Summary")
    st.write(data.get("summary", "No summary available"))

    # Species
    st.subheader("Species")
    species_list = data.get("species", [])
    if species_list:
        for species_group in species_list:
            scientific_name = next(
                (
                    item["value"]
                    for item in species_group
                    if item.get("accession") == "MS:1001469"
                ),
                "Unknown",
            )
            taxid = next(
                (
                    item["value"]
                    for item in species_group
                    if item.get("accession") == "MS:1001467"
                ),
                "N/A",
            )
            st.write(f"**{scientific_name}** (NCBI TaxID: {taxid})")
    else:
        st.write("No species information available")

    # Dataset Links
    st.subheader("Dataset Links")
    dataset_links = data.get("datasetLink", [])
    for link in dataset_links:
        link_name = link.get("name", "Link")
        link_url = link.get("value", "")
        if link_url:
            st.markdown(f"- [{link_name}]({link_url})")

    # Publications
    st.subheader("Publications")
    publications = data.get("publications", [])
    if publications:
        for pub in publications:
            pub_name = pub.get("name", "Unknown publication")
            pub_value = pub.get("value", "")
            if pub_value:
                st.write(f"- {pub_value}")
            else:
                st.write(f"- {pub_name}")
    else:
        st.write("No publications available")


def display_metabolights_metadata(data: Dict[str, Any], dataset_id: str) -> None:
    """Display MetaboLights metadata"""

    # Extract study information
    isa_investigation = data.get("isaInvestigation", {})
    studies = isa_investigation.get("studies", [])

    if not studies:
        st.error("No study information found")
        return

    study = studies[0]  # Get first study
    mtbls_study = data.get("mtblsStudy", {})

    # Title
    st.header(study.get("title", "No title available"))

    # Dataset ID badge
    st.markdown(f"**Dataset ID:** `{dataset_id}`")

    # Study status
    status = mtbls_study.get("studyStatus", "Unknown")
    status_color = "green" if status == "Public" else "orange"
    st.markdown(f"**Status:** :{status_color}[{status}]")

    # Summary/Description
    st.subheader("Summary")
    st.write(study.get("description", "No summary available"))

    # Species/Organisms
    st.subheader("Species")
    organisms = study.get("studyDesignDescriptors", [])
    organism_found = False
    for descriptor in organisms:
        if "organism" in descriptor.get("annotationValue", "").lower():
            st.write(f"**{descriptor.get('annotationValue', 'Unknown')}**")
            organism_found = True

    if not organism_found:
        # Try to extract from study factors or characteristics
        factors = study.get("factors", [])
        for factor in factors:
            if "organism" in factor.get("factorName", "").lower():
                st.write(f"**{factor.get('factorName', 'Unknown')}**")
                organism_found = True
                break

    if not organism_found:
        st.write("Species information not explicitly listed")

    # Dataset Links
    st.subheader("Dataset Links")
    st.markdown(
        f"- [MetaboLights Study Page](https://www.ebi.ac.uk/metabolights/{dataset_id})"
    )
    if mtbls_study.get("studyHttpUrl"):
        st.markdown(f"- [FTP HTTP Access]({mtbls_study.get('studyHttpUrl')})")
    if mtbls_study.get("studyFtpUrl"):
        st.markdown(f"- [FTP Access]({mtbls_study.get('studyFtpUrl')})")

    # Publications
    st.subheader("Publications")
    publications = study.get("studyPublications", [])
    if publications:
        for pub in publications:
            title = pub.get("title", "")
            doi = pub.get("doi", "")
            pubmed_id = pub.get("pubMedID", "")

            if doi:
                st.markdown(f"- [{title}](https://doi.org/{doi})")
            elif pubmed_id:
                st.markdown(
                    f"- [{title}](https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/)"
                )
            elif title:
                st.write(f"- {title}")
    else:
        st.write("No publications available")

    # Additional Information in Expander
    with st.expander("Additional Information"):
        # Submission and Release Dates
        st.write(f"**Submission Date:** {study.get('submissionDate', 'N/A')}")
        st.write(f"**Public Release Date:** {study.get('publicReleaseDate', 'N/A')}")

        # People/Contacts
        people = study.get("people", [])
        if people:
            st.write("**Contacts:**")
            for person in people:
                name = f"{person.get('firstName', '')} {person.get('lastName', '')}".strip()
                email = person.get("email", "")
                affiliation = person.get("affiliation", "")
                roles = person.get("roles", [])
                role_names = [role.get("annotationValue", "") for role in roles]

                if name:
                    role_str = f" ({', '.join(role_names)})" if role_names else ""
                    st.write(f"- **{name}**{role_str}")
                    if affiliation:
                        st.write(f"  - {affiliation}")
                    if email:
                        st.write(f"  - {email}")

        # Assays
        assays = study.get("assays", [])
        if assays:
            st.write("**Assays:**")
            for assay in assays:
                measurement = assay.get("measurementType", {}).get(
                    "annotationValue", "Unknown"
                )
                technology = assay.get("technologyType", {}).get(
                    "annotationValue", "Unknown"
                )
                platform = assay.get("technologyPlatform", "N/A")
                st.write(f"- {measurement} using {technology}")
                if platform != "N/A":
                    st.write(f"  - Platform: {platform}")


def display_workbench_metadata(data: Dict[str, Any], dataset_id: str) -> None:
    """Display Metabolomics Workbench metadata"""

    # Extract sections
    project = data.get("PROJECT", {})
    study = data.get("STUDY", {})
    subject = data.get("SUBJECT", {})
    collection = data.get("COLLECTION", {})

    # Title
    study_title = study.get("STUDY_TITLE") or project.get(
        "PROJECT_TITLE", "No title available"
    )
    st.header(study_title)

    # Dataset ID badge
    st.markdown(f"**Dataset ID:** `{dataset_id}`")

    # Summary
    st.subheader("Summary")
    summary = study.get("STUDY_SUMMARY") or project.get(
        "PROJECT_SUMMARY", "No summary available"
    )
    st.write(summary)

    # Species
    st.subheader("Species")
    species_name = subject.get("SUBJECT_SPECIES", "Unknown")
    taxonomy_id = subject.get("TAXONOMY_ID", "N/A")
    subject_type = subject.get("SUBJECT_TYPE", "Unknown")
    st.write(f"**{species_name}** ({subject_type})")
    if taxonomy_id != "N/A":
        st.write(f"NCBI TaxID: {taxonomy_id}")

    # Dataset Links
    st.subheader("Dataset Links")
    st.markdown(
        f"- [Metabolomics Workbench Study Page](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID={dataset_id})"
    )

    # Additional Information in Expander
    with st.expander("Additional Information"):
        # Study Information
        st.write("**Study Information:**")
        st.write(
            f"- Institute: {study.get('INSTITUTE', project.get('INSTITUTE', 'N/A'))}"
        )
        st.write(f"- Total Subjects: {study.get('TOTAL_SUBJECTS', 'N/A')}")
        num_males = study.get("NUM_MALES", "N/A")
        num_females = study.get("NUM_FEMALES", "N/A")
        st.write(f"- Males: {num_males}, Females: {num_females}")

        # Contact Information
        st.write("**Contact:**")
        first_name = study.get("FIRST_NAME") or project.get("FIRST_NAME", "")
        last_name = study.get("LAST_NAME") or project.get("LAST_NAME", "")
        email = study.get("EMAIL") or project.get("EMAIL", "")
        if first_name or last_name:
            st.write(f"- {first_name} {last_name}")
        if email:
            st.write(f"- {email}")

        # Subject Information
        if subject:
            st.write("**Subject Information:**")
            age_range = subject.get("AGE_OR_AGE_RANGE", "N/A")
            gender = subject.get("GENDER", "N/A")
            if age_range != "N/A":
                st.write(f"- Age Range: {age_range}")
            if gender != "N/A":
                st.write(f"- Gender: {gender}")

            lifestyle = subject.get("HUMAN_LIFESTYLE_FACTORS", "")
            if lifestyle and lifestyle != "NA":
                st.write(f"- Lifestyle Factors: {lifestyle}")

        # Collection Information
        if collection:
            st.write("**Sample Collection:**")
            sample_type = collection.get("SAMPLE_TYPE", "N/A")
            location = collection.get("COLLECTION_LOCATION", "N/A")
            storage = collection.get("STORAGE_CONDITIONS", "N/A")

            if sample_type != "N/A":
                st.write(f"- Sample Type: {sample_type}")
            if location != "N/A":
                st.write(f"- Collection Location: {location}")
            if storage != "N/A":
                st.write(f"- Storage Conditions: {storage}")


def display_workbench_metadata_summary(data: Dict[str, Any], dataset_id: str) -> None:
    """Display only summary information for workbench metadata."""
    st.subheader("Summary Information")
    st.write(f"**Dataset ID**: `{dataset_id}`")
    st.write(f"**Study Title**: {data.get('study_title', 'N/A')}")
    st.write(f"**Species**: {data.get('species', 'N/A')}")
    st.write(f"**Institute**: {data.get('institute', 'N/A')}")
    st.write(f"**Analysis Type**: {data.get('analysis_type', 'N/A')}")
    st.write(f"**Number of Samples**: {data.get('number_of_samples', 'N/A')}")
    st.write(f"**[Metabolomics Workbench Link](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?StudyID={dataset_id})**")


def display_metabolomics_metadata(dataset_id: str) -> None:
    """
    Main function to fetch and display metadata for a metabolomics dataset.
    Automatically detects the repository based on dataset ID prefix.

    Parameters:
    -----------
    dataset_id : str
        The dataset identifier (e.g., 'MSV000082221', 'MTBLS916', 'ST002320')
    """

    dataset_id = dataset_id.strip()

    # Detect repository
    repository = detect_repository(dataset_id)

    if repository is None:
        st.error(f"Unrecognized dataset ID format: {dataset_id}")
        st.info(
            "Supported formats: MSV* (MassIVE), MTBLS* (MetaboLights), ST* (Metabolomics Workbench)"
        )
        return

    try:
        # Show loading indicator
        with st.spinner(f"Fetching metadata from {repository.value}..."):

            # Fetch data based on repository
            if repository == DataRepository.MASSIVE:
                data = fetch_massive_metadata(dataset_id)
                display_massive_metadata(data, dataset_id)

            elif repository == DataRepository.METABOLIGHTS:
                data = fetch_metabolights_metadata(dataset_id)
                st.write(data)
                display_metabolights_metadata(data, dataset_id)

            elif repository == DataRepository.WORKBENCH:
                try:
                    data = fetch_workbench_metadata(dataset_id)
                    display_workbench_metadata(data, dataset_id)
                except requests.JSONDecodeError as e:
                    data = fetch_workbench_metadata(dataset_id, summary_only=True)
                    display_workbench_metadata_summary(data, dataset_id)

        # Add repository badge
        st.markdown(f"---")
        st.caption(f"Data source: {repository.value.upper()}")

    except requests.exceptions.RequestException as e:
        st.error(f"Error fetching data for {dataset_id}: {str(e)}")
    except KeyError as e:
        st.error(f"Error parsing data - missing expected field: {str(e)}")
    except Exception as e:
        raise
        st.error(f"Error processing data: {str(e)}")
