import streamlit as st
import pandas as pd
import requests


FALLBACK_OPTIONS = [
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


@st.cache_data(ttl=3600)
def fetch_database_options():
    url = "https://fasst.gnps2.org/libraries"
    response = requests.get(url, timeout=5)
    response.raise_for_status()
    libraries = response.json()
    return [lib['value'] for lib in libraries]


def get_database_options():
    with st.spinner('Retrieving MASST Database options...'):
        try:
            return fetch_database_options()
        except requests.Timeout:
            st.toast("Database Request timed out. Using fallback database options.", icon="⚠️")
        except requests.RequestException as e:
            st.toast(f"Error fetching database options: {e}", icon="❌")
            st.toast("Using fallback database options.", icon="⚠️")
        return FALLBACK_OPTIONS


def create_masst_sidebar():
    """Create sidebar widgets for MASST parameters and return the values"""

    with st.sidebar:
        with st.expander("MASST Query Parameters", icon=':material/settings:', expanded=False):

            # Database selection
            options_list = sorted(get_database_options())
            database = st.selectbox(
                "Database",
                options=options_list,
                index=options_list.index("metabolomicspanrepo_index_latest") if "metabolomicspanrepo_index_latest" in options_list else 0,
                help="Type of database to search"
            )

            # # MASST type
            # masst_type = st.selectbox(
            #     "MASST Type",
            #     options=["masst", "gnpsdata", "microbemasst"],
            #     index=0,
            #     help="Type of MASST to give results"
            # )

            # Analog search
            analog = st.selectbox(
                "Analog Search",
                options=["No", "Yes"],
                index=0,
                help="Perform analog search"
            )

            # Tolerance parameters
            st.subheader("Tolerance Parameters")

            precursor_tolerance = st.number_input(
                "Precursor m/z Tolerance",
                min_value=0.001,
                max_value=1.0,
                value=0.02,
                step=0.001,
                format="%.3f",
                help="Precursor mass tolerance"
            )

            fragment_tolerance = st.number_input(
                "Fragment m/z Tolerance",
                min_value=0.001,
                max_value=1.0,
                value=0.02,
                step=0.001,
                format="%.3f",
                help="Fragment mass tolerance"
            )

            # Cosine threshold
            cosine_threshold = st.number_input(
                "Cosine Similarity Threshold",
                min_value=0.0,
                max_value=1.0,
                value=0.7,
                step=0.01,
                help="Minimum cosine similarity score"
            )

        return {
            'database': database,
            'masst_type': "masst",
            'analog': analog == "Yes",
            'precursor_mz_tol': precursor_tolerance,
            'fragment_mz_tol': fragment_tolerance,
            'min_cos': cosine_threshold
        }


def create_usi_input():
    """Create main area input for USI data"""
    with st.sidebar:
        st.subheader("USI Input Data", help='Enter one USI per line in the data editor below. You can add more rows as necessary.' )

        # Create sample data structure
        if 'usi_data' not in st.session_state:
            st.session_state.usi_data = pd.DataFrame({
                'usi': [
                    'USI-1',
                    'USI-2',
                    'USI-3'
                ],
                'compound_name': ['Compound-name-1', 'Compound-name-2', 'Compound-name-3']
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
                    width="small"
                ),
                "compound_name": st.column_config.TextColumn(
                    "Compound Name",
                    help="Name or identifier for the compound",
                    width="medium"
                )
            }
        )

    return edited_data