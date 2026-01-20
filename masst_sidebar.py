import streamlit as st
import pandas as pd
import requests

from utils.query_params import sync_param, get_default


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


def create_masst_sidebar(disabled: bool = False, initial_params: dict = None) -> dict:
    """Create sidebar widgets for MASST parameters and return the values.
    
    Args:
        disabled: Whether to disable all widgets
        initial_params: Optional dict with initial values from query params
                       Keys: database, analog, precursor_tol, fragment_tol, min_cos
    """
    if initial_params is None:
        initial_params = {}

    with st.sidebar:
        with st.expander("MASST Query Parameters", icon=':material/settings:', expanded=False):

            # Database selection
            options_list = sorted(get_database_options())
            
            # Get initial database value
            initial_db = initial_params.get('database', get_default('database'))
            if initial_db not in options_list:
                initial_db = get_default('database')
            db_index = options_list.index(initial_db) if initial_db in options_list else 0
            
            database = st.selectbox(
                "Database",
                options=options_list,
                index=db_index,
                help="Type of database to search",
                disabled=disabled,
                key='masst_database'
            )

            # # MASST type
            # masst_type = st.selectbox(
            #     "MASST Type",
            #     options=["masst", "gnpsdata", "microbemasst"],
            #     index=0,
            #     help="Type of MASST to give results"
            # )

            # Analog search
            initial_analog = initial_params.get('analog', get_default('analog'))
            analog = st.selectbox(
                "Analog Search",
                options=["No", "Yes"],
                index=1 if initial_analog else 0,
                help="Perform analog search",
                disabled=disabled,
                key='masst_analog'
            )

            # Tolerance parameters
            st.subheader("Tolerance Parameters")

            initial_precursor = initial_params.get('precursor_tol', get_default('precursor_tol'))
            precursor_tolerance = st.number_input(
                "Precursor m/z Tolerance",
                min_value=0.001,
                max_value=1.0,
                value=float(initial_precursor),
                step=0.001,
                format="%.3f",
                help="Precursor mass tolerance",
                disabled=disabled,
                key='masst_precursor_tol'
            )

            initial_fragment = initial_params.get('fragment_tol', get_default('fragment_tol'))
            fragment_tolerance = st.number_input(
                "Fragment m/z Tolerance",
                min_value=0.001,
                max_value=1.0,
                value=float(initial_fragment),
                step=0.001,
                format="%.3f",
                help="Fragment mass tolerance",
                disabled=disabled,
                key='masst_fragment_tol'
            )

            # Cosine threshold
            initial_cos = initial_params.get('min_cos', get_default('min_cos'))
            cosine_threshold = st.number_input(
                "Cosine Similarity Threshold",
                min_value=0.0,
                max_value=1.0,
                value=float(initial_cos),
                step=0.01,
                help="Minimum cosine similarity score",
                disabled=disabled,
                key='masst_min_cos'
            )

        # Sync parameters to URL in real-time
        sync_param('database', database)
        sync_param('analog', analog == "Yes")
        sync_param('precursor_tol', precursor_tolerance)
        sync_param('fragment_tol', fragment_tolerance)
        sync_param('min_cos', cosine_threshold)

        return {
            'database': database,
            'masst_type': "masst",
            'analog': analog == "Yes",
            'precursor_mz_tol': precursor_tolerance,
            'fragment_mz_tol': fragment_tolerance,
            'min_cos': cosine_threshold
        }


def create_usi_input(disabled: bool = False, usi_data: pd.DataFrame = None, sync_url: bool = True) -> pd.DataFrame:
    """Create main area input for USI data.
    
    Args:
        disabled: Whether to disable the widget
        usi_data: Initial DataFrame with 'usi' and 'compound_name' columns
        sync_url: Whether to sync changes to URL query params
    """
    with st.sidebar:
        st.subheader("USI Input Data", help='Enter one USI per line in the data editor below. You can add more rows as necessary.' )

        # Create sample data structure
        if not isinstance(usi_data, pd.DataFrame):
            usi_data = pd.DataFrame({
                'usi': [
                    'USI-1',
                    'USI-2',
                    'USI-3'
                ],
                'compound_name': ['Compound-name-1', 'Compound-name-2', 'Compound-name-3']
            })

        # Data editor
        edited_data = st.data_editor(
            usi_data,
            disabled=disabled,
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
            },
            key='usi_data_editor'
        )
        
        # Sync USI data to URL (only if not using example data)
        if sync_url and not st.session_state.get('use_example', False):
            sync_param('usi', edited_data)

    return edited_data