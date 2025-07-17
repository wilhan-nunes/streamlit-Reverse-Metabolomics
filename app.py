import os

import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import io
import subprocess

from welcome import welcome_page


def get_git_short_rev():
    try:
        return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], stderr=subprocess.DEVNULL).decode().strip()
    except subprocess.CalledProcessError:
        return ".git/ not found"

# Configure matplotlib for PDF output
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# TODO: Bump version
app_version = "2025-07-17"
git_hash = get_git_short_rev()
repo_link = "https://github.com/wilhan-nunes/streamlit-Reverse-Metabolomics"

st.set_page_config(
    page_title="Reverse Metabolomics Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={"About": (f"**App version**: {app_version} | "
                          f"[**Git Hash**: {git_hash}]({repo_link}/commit/{git_hash})")},
)

st.title("🧬 Reverse Metabolomics Analysis Tool")
st.markdown("Edit entries in the sidebar to analyze fastMASST results with interactive visualizations")


def check_memory_usage():
    """Check current memory usage and warn if high"""
    import psutil
    memory_percent = psutil.virtual_memory().percent
    if memory_percent > 80:
        st.warning(f"High memory usage detected: {memory_percent:.1f}%. Consider reducing dataset size.")
    return memory_percent


# Create custom colormap
@st.cache_data
def create_colormap():
    colors = [(1, 1, 1), (0.78, 0.84, 0.94), (0.92, 0.69, 0.65)]  # white -> blue -> red
    n_bins = 100
    cmap_name = 'white_blue_red'
    return LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)


cmap_wbr = create_colormap()


# Helper functions
def analyze_counts(df, column_interest):
    """Prepare a table with counts of the fastMASST results."""
    list_body_parts = df[column_interest].unique().tolist()
    df_body_parts = pd.DataFrame(list_body_parts, columns=[column_interest])

    df_counts = df[column_interest].value_counts().rename_axis(column_interest).reset_index(name='Counts_fastMASST')

    compounds = df.groupby(column_interest)['Compound'].agg(
        ['nunique', lambda x: ', '.join(map(str, x.unique()))]).reset_index()
    compounds.columns = [column_interest, 'Compounds', 'CompoundsList']

    combined = pd.merge(df_body_parts, df_counts, on=column_interest, how='left')
    combined = pd.merge(combined, compounds, on=column_interest, how='left')

    return combined


def prepare_pivot_table(df, column_interest, compound, log_transform=False, normalize_redu=False,
                        df_redu_counts=None, redu_count_col=None):
    """Prepare pivot table for heatmap visualization."""
    grouped_df = df.groupby([compound, column_interest]).size().reset_index(name='Count')
    pivot_table = grouped_df.pivot(index=column_interest, columns=compound, values='Count').fillna(0)
    pivot_table.reset_index(inplace=True)

    if normalize_redu and df_redu_counts is not None:
        # Merge with ReDU counts for normalization
        pivot_table = pd.merge(pivot_table, df_redu_counts, on=column_interest, how='left')

        # Normalize by ReDU counts
        columns_to_normalize = pivot_table.columns.difference([column_interest, redu_count_col])
        pivot_table[columns_to_normalize] = pivot_table[columns_to_normalize].div(pivot_table[redu_count_col], axis=0)
        pivot_table.drop(redu_count_col, axis=1, inplace=True)

        # Calculate relative abundance
        sums = pivot_table.select_dtypes(include='number').sum().to_dict()
        sums[column_interest] = 'Sum'
        pivot_table = pd.concat([pivot_table, pd.DataFrame([sums])], ignore_index=True)

        for column in pivot_table.select_dtypes(include='number').columns:
            total_value = pivot_table[column].iloc[-1]
            if total_value > 0:
                pivot_table[column] = pivot_table[column] / total_value * 100

        pivot_table = pivot_table[pivot_table[column_interest] != 'Sum']

    if log_transform:
        columns_for_transform = pivot_table.columns[1:]
        pivot_table[columns_for_transform] = pivot_table[columns_for_transform] + 1
        log_transformed = np.log2(pivot_table[columns_for_transform])
        pivot_table = pd.concat([pivot_table[column_interest], log_transformed], axis=1)

    return pivot_table


@st.cache_data
def load_and_process_data(fastmasst_results: pd.DataFrame, usis_table: pd.DataFrame, df_redu: pd.DataFrame,
                          tolerance: float, organism_filter: str = "All organisms"):
    """Load and process data with pre-loaded ReDU dataframe"""

    usi_to_name = dict(zip(usis_table['usi'], usis_table['compound_name']))
    compound_names = list(usis_table['compound_name'])

    df_combined = fastmasst_results.copy()
    # filter for tolerance
    df_combined = df_combined[(df_combined['Delta Mass'] >= -tolerance) & (df_combined['Delta Mass'] <= tolerance)]
    df_combined['Compound'] = df_combined['query_usi'].apply(lambda x: usi_to_name.get(x, 'Unknown'))

    # Create filepath column
    df_combined['filepath'] = df_combined['Dataset'] + "/" + df_combined['USI'].str.split('/').str[-1]
    df_combined['filepath'] = df_combined['filepath'].str.split(':').str[0]
    df_combined['filepath'] = df_combined['filepath'].str.replace('.mzML', '').str.replace('.mzXML', '')

    # Merge datasets (df_redu is already processed)
    df_merged = pd.merge(df_combined, df_redu, left_on='filepath', right_on='filepath', how='left',
                         suffixes=('_fasst', '_redu'))

    # Standardize body part names
    if 'UBERONBodyPartName' in df_merged.columns:
        body_part_replacements = {
            'skin of trunk': 'skin',
            'skin of pes': 'skin',
            'head or neck skin': 'skin',
            'axilla skin': 'skin',
            'skin of manus': 'skin',
            'arm skin': 'skin',
            'blood plasma': 'blood',
            'blood serum': 'blood'
        }

        for old, new in body_part_replacements.items():
            df_merged['UBERONBodyPartName'] = df_merged['UBERONBodyPartName'].str.replace(old, new)

    taxonomy_filters = {
        "Humans": ['9606|Homo sapiens'],
        "Rodents": ['10088|Mus', '10090|Mus musculus', '10105|Mus minutoides', '10114|Rattus',
                    '10116|Rattus norvegicus']
    }

    if organism_filter in taxonomy_filters:
        selected_taxa = taxonomy_filters[organism_choice]
        df_filtered = df_merged[df_merged['NCBITaxonomy'].isin(selected_taxa)]
        df_redu_filtered = df_redu[df_redu['NCBITaxonomy'].isin(selected_taxa)]
        return df_filtered, df_redu_filtered

    else:
        return df_merged, df_redu



def create_heatmap(pivot_table, variable, metric, log_scale=False, normalize=False,
                   col_cluster=False, row_cluster=False, width=2, height=8, ):
    """Create heatmap visualization."""
    columns_for_heatmap = pivot_table.columns[1:]

    fig = sns.clustermap(
        pivot_table[columns_for_heatmap],
        col_cluster=col_cluster,
        row_cluster=row_cluster,
        method='complete',
        metric=metric,
        cmap=cmap_wbr,
        yticklabels=pivot_table[variable],
        linewidths=0.005,
        linecolor='white',
        cbar_kws={'orientation': 'horizontal'},
        figsize=(width, height)
    )

    fig.ax_heatmap.yaxis.set_ticks_position('left')
    fig.ax_heatmap.yaxis.set_label_position('left')
    fig.ax_heatmap.set_xticklabels(fig.ax_heatmap.get_xticklabels(), rotation=90)

    cax = fig.ax_cbar
    cax.set_position((-0.25, -0.15, .3, .05))
    cbar = fig.ax_heatmap.collections[0].colorbar

    if log_scale:
        cbar.set_label('log2(spectral matches)', fontsize=10)
    elif normalize:
        cbar.set_label('Normalized spectral matches (%)', fontsize=10)
    else:
        cbar.set_label('Spectral matches', fontsize=10)

    fig.ax_heatmap.xaxis.set_label_coords(0.5, -0.5)
    heatmap_title = f"{variable}"
    fig.ax_heatmap.set_title(heatmap_title, fontsize=14, weight='bold')

    line_count = len(pivot_table[variable].unique())
    fig.ax_heatmap.axhline(y=line_count, color='black', linewidth=1.5)
    fig.ax_heatmap.axvline(x=0, color='black', linewidth=1.5)

    return fig

# Cache the ReDU dataframe loading
@st.cache_data
def load_redu_data():
    """Load ReDU metadata file and return as DataFrame"""
    try:
        df_redu = pd.read_csv('REDU_metadata.tsv', sep='\t')
        # Process ReDU table
        df_redu['filename_2'] = df_redu['filename'].str.split('/').str[-1]
        df_redu['filename_2'] = df_redu['filename_2'].str.replace('.mzML', '').str.replace('.mzXML', '')
        df_redu['filepath'] = df_redu['ATTRIBUTE_DatasetAccession'].astype(str) + '/' + df_redu['filename_2'].astype(str)
        return df_redu
    except Exception as e:
        st.error(f"Error loading ReDU metadata: {str(e)}")
        return None



# Sidebar for inputs
with st.sidebar:
    st.header("📁 Data Input")

    # Example data option
    use_example = st.checkbox("Load example data", value=False, help="Use built-in example files instead of uploading your own", key='use_example')

    if not use_example:
        from masst_sidebar import create_masst_sidebar, create_usi_input
        from masst_client import masst_query_all

        usi_data = create_usi_input()
        masst_query_params = create_masst_sidebar()
        # Filter out empty rows
        query_df = usi_data[usi_data['usi'].str.strip() != ''].copy()

        st.session_state.masst_query_params = masst_query_params
        st.session_state['run_masst_query'] = st.button("Run MASST Query", type="primary")

    else:
        results = pd.read_csv('example_data/his-c4-phe-c4-phe-ca.csv')
        usi_data = pd.DataFrame({
            'usi': [
                'mzspec:gnps:GNPS-LIBRARY:accession:CCMSLIB00006582001',
                'mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00010010601',
                'mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00011434738'
            ],
            'compound_name': ['Phe-CA', 'Phe-C4:0', 'His-C4:0']
        })
        st.write(usi_data)
        mass_tolerance = st.number_input("Delta mass tolerance (Da)", min_value=0.0, max_value=1.0,
                                                 value=0.02, step=0.01)
        if st.button('Run Example'):
            st.session_state.results = results

    if st.button("Clear All Cache", type="secondary"):
        st.session_state.clear()
        st.session_state["use_example"] = False
        st.rerun()

# Download ReDU metadata if not exists
if not os.path.exists('REDU_metadata.tsv'):
    from download_redu import download_redu_metadata
    with st.spinner("Downloading ReDU metadata file... this may take a while!"):
        download_redu_metadata('REDU_metadata.tsv')




if st.session_state.get('run_masst_query', False):
    if len(query_df) > 0:
        with st.spinner("Running MASST query..."):
            # Here you would call your imported masst_query_all function
            results = masst_query_all(query_df, **masst_query_params)
            st.success(f"Query performed with {len(query_df)} USIs with parameters: {masst_query_params}")
            st.session_state.results = results
            # st.write(results)
    else:
        st.error("Please add at least one USI to query")

if "results" in st.session_state:
    # Load ReDU data once and cache it
    df_redu = load_redu_data()
    if df_redu is None:
        st.error("Failed to load ReDU metadata file")
        st.stop()

    results = st.session_state.results

    try:
        if use_example:
            st.success("Using example data filtered for Precursor delta mass tolerance of %.2f Da" % mass_tolerance)
        else:
            masst_query_params = st.session_state.masst_query_params
            mass_tolerance = masst_query_params.get('precursor_mz_tol', 0.02)

               # Organism selection
        col1, col2 = st.columns(2)
        with col1:
            st.write("🧬 Organism Filter")

            # available_organisms = df_merged['NCBITaxonomy'].dropna().unique()

            organism_choice = st.selectbox(
                "Select organism:",
                options=['Humans', 'Rodents', 'All organisms'],
                help="Filter data by organism type"
            )

            df_filtered, df_redu_filtered = load_and_process_data(results, usi_data, df_redu, mass_tolerance, organism_choice)

        with col2:
            # Analysis options
            st.write("📈 Analysis Options")

            available_columns = [col for col in df_filtered.columns if
                                 col in ['UBERONBodyPartName', 'DOIDCommonName', 'ATTRIBUTE_SubjectGender',
                                         'ATTRIBUTE_Age']]

            analysis_column = st.selectbox(
                "Select analysis variable:",
                options=available_columns,
                help="Choose which ReDU column to analyze"
            )

            unique_options = list(df_filtered[analysis_column].unique())
            variable_to_exclude = st.multiselect("Exclude from analysis (optional)", options=unique_options)

        st.info(f"📊 Filtered data: {len(df_filtered):,} spectral matches for {organism_choice.lower()}")

        # Main content
        col1, col2 = st.columns([1, 2])

        with col1:
            st.header("📋 Summary Statistics")

            # Counts table
            if analysis_column:
                counts_table = analyze_counts(df_filtered, analysis_column)
                st.subheader(f"Counts by {analysis_column}")
                st.dataframe(counts_table)

                # Download counts table
                csv_counts = counts_table.to_csv(index=False)
                st.download_button(
                    label="📥 Download counts table",
                    data=csv_counts,
                    file_name=f"counts_{analysis_column}.csv",
                    mime="text/csv"
                )

        with col2:
            st.header("🔥 Heatmap Visualizations")

            # Heatmap options
            col1, col2 = st.columns(2)
            with col1:
                heatmap_type = st.selectbox(
                    "Heatmap type:",
                    options=['Raw counts', 'Log-transformed counts', 'ReDU-normalized counts'],
                    help="Choose the type of heatmap to generate"
                )

            clustering_metrics = {
                "Bray-Curtis": "braycurtis",
                "Canberra": "canberra",
                "Chebyshev": "chebyshev",
                "City Block (Manhattan)": "cityblock",
                "Correlation": "correlation",
                "Cosine": "cosine",
                "Dice": "dice",
                "Euclidean": "euclidean",
                "Hamming": "hamming",
                "Jaccard": "jaccard",
                "Jensen-Shannon": "jensenshannon",
                "Kulczynski": "kulczynski1",
                "Mahalanobis": "mahalanobis",
                "Matching": "matching",
                "Minkowski": "minkowski",
                "Rogers-Tanimoto": "rogerstanimoto",
                "Russell-Rao": "russellrao",
                "Standardized Euclidean": "seuclidean",
                "Sokal-Michener": "sokalmichener",
                "Sokal-Sneath": "sokalsneath",
                "Squared Euclidean": "sqeuclidean",
                "Yule": "yule"
            }

            with col2:
                metric = st.selectbox(
                    "Clustering metric:",
                    options=list(clustering_metrics.keys()),
                    help='The distance metric to use. See scipy.spatial.distance.pdist documentation for more info.'
                )

            col1, col2 = st.columns(2)

            with col1:
                col_cluster = st.checkbox("Cluster columns", value=False)
                heatmap_width = st.number_input("Heatmap width", min_value=2, max_value=10, value=4)
            with col2:
                row_cluster = st.checkbox("Cluster rows", value=False)
                heatmap_height = st.number_input("Heatmap height", min_value=4, max_value=15, value=4)

            if st.button("🎨 Generate Heatmap"):
                with st.spinner("Generating heatmap..."):
                    try:
                        log_transform = heatmap_type == 'Log-transformed counts'
                        normalize_redu = heatmap_type == 'ReDU-normalized counts'

                        df_redu_counts = None
                        redu_count_col = None

                        if normalize_redu:
                            df_redu_counts = df_redu_filtered[analysis_column].value_counts().reset_index()
                            df_redu_counts.columns = [analysis_column, f'{analysis_column}_counts']
                            redu_count_col = f'{analysis_column}_counts'

                        if variable_to_exclude:
                            drop_mask = df_filtered[analysis_column].isin(variable_to_exclude)
                            indices_to_drop = df_filtered[drop_mask].index
                            df_filtered = df_filtered.drop(indices_to_drop)

                        pivot_table = prepare_pivot_table(
                            df_filtered,
                            analysis_column,
                            'Compound',
                            log_transform=log_transform,
                            normalize_redu=normalize_redu,
                            df_redu_counts=df_redu_counts,
                            redu_count_col=redu_count_col
                        )

                        fig = create_heatmap(
                            pivot_table,
                            analysis_column,
                            metric=clustering_metrics.get(metric, 'euclidean'),
                            log_scale=log_transform,
                            normalize=normalize_redu,
                            col_cluster=col_cluster,
                            row_cluster=row_cluster,
                            width=heatmap_width,
                            height=heatmap_height
                        )

                        st.pyplot(fig.figure)

                        # Download options
                        col_a, col_b = st.columns(2)

                        with col_a:
                            # Download pivot table
                            csv_pivot = pivot_table.to_csv(index=False)
                            st.download_button(
                                label="📥 Download pivot table",
                                data=csv_pivot,
                                file_name=f"pivot_table_{analysis_column}_{heatmap_type.replace(' ', '_').lower()}.csv",
                                mime="text/csv"
                            )

                        with col_b:
                            # Download plot as PNG
                            img_buffer = io.BytesIO()
                            fig.figure.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
                            img_buffer.seek(0)

                            st.download_button(
                                label="📥 Download heatmap (PNG)",
                                data=img_buffer.getvalue(),
                                file_name=f"heatmap_{analysis_column}_{heatmap_type.replace(' ', '_').lower()}.png",
                                mime="image/png"
                            )

                    except Exception as e:
                        st.error(f"Error generating heatmap: {str(e)}")

        # Data preview
        st.header("Data Preview")
        st.subheader("Merged Dataset", help="Showing first 100 rows of the merged dataset")
        st.dataframe(df_filtered.head(100))

        # Download full dataset
        csv_full = df_filtered.to_csv(index=False)
        st.download_button(
            label="📥 Download full merged dataset",
            data=csv_full,
            file_name=f"merged_dataset_{organism_choice.lower().replace(' ', '_')}.csv",
            mime="text/csv"
        )

    except Exception as e:
        st.error(f"❌ Error processing data: {str(e)}")
        st.info("Please check your file formats and ensure they match the expected structure.")
        raise

else:
    st.info("👈 Please edit the USI input data table in the sidebar, choose your parameters and click Run MASST Query to start.")

    welcome_page()