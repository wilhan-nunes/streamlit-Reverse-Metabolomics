import streamlit as st
import pandas as pd
from masst_client import masst_query_all  # Assuming masst_client.py is in the same directory
import requests


@st.cache_data
def load_redu_data():
    """Load ReDU metadata - you'll need to provide the actual path or URL"""
    # This is a placeholder - replace with actual ReDU data loading
    # You might want to download it from a URL or load from a file
    try:
        # Example: Load from a URL or local file
        df_redu = pd.read_csv("example_data/REDU_metadata_nat_prot.tsv", sep='\t')

        # For now, return None to indicate ReDU data is not available
        return df_redu
    except:
        return None


def get_database_options():
    url = "https://fasst.gnps2.org/libraries"
    try:
        response = requests.get(url)
        response.raise_for_status()
        libraries = response.json()
        return [lib['value'] for lib in libraries]
    except requests.RequestException as e:
        st.error(f"Error fetching database options: {e}")
        return [
            "gnpsdata_index",
            "gnpslibrary",
            "massivedata_index",
            "massivekb_index",
            "metabolomicspanrepo_index_latest",
            "metabolomicspanrepo_index_nightly",
            "panrepo_2024_11_12",
            "panrepo_2025_06_18",
            "ptfi2_index",
            "gnpsdata_test_index",
            "ORNL_Bioscales2",
            "ORNL_Populus_LC_MSMS"
        ]


def create_masst_sidebar():
    """Create sidebar widgets for MASST parameters and return the values"""

    st.sidebar.header("MASST Query Parameters")

    # Database selection
    options_list = sorted(get_database_options())
    database = st.sidebar.selectbox(
        "Database",
        options=options_list,
        index=options_list.index("panrepo_2024_11_12") if "panrepo_2024_11_12" in options_list else 0,
        help="Type of database to search"
    )

    # MASST type
    masst_type = st.sidebar.selectbox(
        "MASST Type",
        options=["masst", "gnpsdata", "microbemasst"],
        index=0,
        help="Type of MASST to give results"
    )

    # Analog search
    analog = st.sidebar.selectbox(
        "Analog Search",
        options=["No", "Yes"],
        index=0,
        help="Perform analog search"
    )

    # Tolerance parameters
    st.sidebar.subheader("Tolerance Parameters")

    precursor_tolerance = st.sidebar.number_input(
        "Precursor m/z Tolerance",
        min_value=0.001,
        max_value=1.0,
        value=0.02,
        step=0.001,
        format="%.3f",
        help="Precursor mass tolerance"
    )

    fragment_tolerance = st.sidebar.number_input(
        "Fragment m/z Tolerance",
        min_value=0.001,
        max_value=1.0,
        value=0.02,
        step=0.001,
        format="%.3f",
        help="Fragment mass tolerance"
    )

    # Cosine threshold
    cosine_threshold = st.sidebar.slider(
        "Cosine Similarity Threshold",
        min_value=0.0,
        max_value=1.0,
        value=0.7,
        step=0.01,
        help="Minimum cosine similarity score"
    )

    return {
        'database': database,
        'masst_type': masst_type,
        'analog': analog == "Yes",
        'precursor_mz_tol': precursor_tolerance,
        'fragment_mz_tol': fragment_tolerance,
        'min_cos': cosine_threshold
    }


def create_usi_input():
    """Create main area input for USI data"""
    with st.sidebar:
        st.subheader("USI Input Data")

        # Create sample data structure
        if 'usi_data' not in st.session_state:
            st.session_state.usi_data = pd.DataFrame({
                'usi': [
                    'mzspec:gnps:GNPS-LIBRARY:accession:CCMSLIB00006582001',
                    'mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00010010601',
                    'mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00011434738'
                ],
                'compound_name': ['Phe-CA', 'Phe-C4:0', 'His-C4:0']
            })

        # Data editor
        edited_data = st.data_editor(
            st.session_state.usi_data,
            num_rows="dynamic",
            use_container_width=True,
            column_config={
                "usi": st.column_config.TextColumn(
                    "USI",
                    help="Universal Spectrum Identifier",
                    width="medium"
                ),
                "compound_name": st.column_config.TextColumn(
                    "Compound Name",
                    help="Name or identifier for the compound",
                    width="medium"
                )
            }
        )

    return edited_data