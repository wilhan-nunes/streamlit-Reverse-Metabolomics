import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.spatial.distance as ssd
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import io
from streamlit.components.v1 import html
from welcome import welcome_page
from download_redu import load_redu
import time
from matplotlib import pyplot as plt

from masst_sidebar import create_masst_sidebar, create_usi_input
from masst_client import masst_query_all
from utils.query_params import load_all_params, sync_param, get_default, clear_all_params, get_param_schema
from utils.retrieve_dataset_metadata import display_metabolomics_metadata


def get_git_short_rev():
    try:
        with open('.git/logs/HEAD', 'r') as f:
            last_line = f.readlines()[-1]
            hash_val = last_line.split()[1]
        return hash_val[:7]
    except Exception:
        return ".git/ not found"

@st.cache_data(ttl=2592000)  # Cache for 30 days
def load_redu_data():
    """Load ReDU metadata file and return as DataFrame"""
    try:
        df, file_date = load_redu(max_age_days=30) # if the file is older than 30 days it will be re-downloaded
        return df, file_date
    except Exception as e:
        st.error(f"Error loading ReDU metadata: {str(e)}")
        return None, None

# Configure matplotlib for PDF output
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# TODO: Bump version
app_version = "2025-10-29"
git_hash = get_git_short_rev()
repo_link = "https://github.com/wilhan-nunes/streamlit-Reverse-Metabolomics"

# Add a tracking token
html('<script async defer data-website-id="74bc9983-13c4-4da0-89ae-b78209c13aaf" src="https://analytics.gnps2.org/umami.js"></script>', width=0, height=0)
html('<script defer src="https://analytics-api.gnps2.org/script.js" data-website-id="74665d88-3b9d-4812-b8fc-7f55ceb08f11"></script>', width=0, height=0)


# Original paper example - https://doi.org/10.1038/s41596-024-01136-2
EXAMPLE_SET_1 = {
                'usi': [
                    'mzspec:gnps:GNPS-LIBRARY:accession:CCMSLIB00006582001',
                    'mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00010010601',
                    'mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00011434738'
                ],
                'compound_name': ['Phe-CA', 'Phe-C4:0', 'His-C4:0']
            }

# HIV related examples
EXAMPLE_SET_2 = {
                'usi': [
                    'mzspec:GNPS2:TASK-fa064fe728814f439a1cd3b72deffcd0-nf_output/clustering/spectra_reformatted.mgf:scan:341',
                    'mzspec:GNPS2:TASK-fa064fe728814f439a1cd3b72deffcd0-nf_output/clustering/spectra_reformatted.mgf:scan:1984',
                    'mzspec:GNPS2:TASK-fa064fe728814f439a1cd3b72deffcd0-nf_output/clustering/spectra_reformatted.mgf:scan:9826',
                    'mzspec:GNPS2:TASK-fa064fe728814f439a1cd3b72deffcd0-nf_output/clustering/spectra_reformatted.mgf:scan:839'
                ],
                'compound_name': ['Histamine-C2:0', 'N-acetylcadaverine-C2:0', 'Candidate-1', 'Candidate-2']
            }


st.set_page_config(
    page_title="Reverse Metabolomics Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={"About": (f"**App version**: {app_version} | "
                          f"[**Git Hash**: {git_hash}]({repo_link}/commit/{git_hash})")},
)


# Load ReDU data directly from cache (shared across all users)
with st.spinner("Loading ReDU metadata..."):
    df_redu, file_date = load_redu_data()

st.title("🧬 Reverse Metabolomics Analysis Tool")


# Create custom colormap
@st.cache_data
def create_colormap():
    colors = [(1, 1, 1), (0.78, 0.84, 0.94), (0.92, 0.69, 0.65)]  # white -> blue -> red
    n_bins = 100
    cmap_name = 'white_blue_red'
    return LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

cmap_wbr = create_colormap()

@st.cache_data
def cached_masst_results(query_df, database, masst_type, analog, precursor_mz_tol, fragment_mz_tol, min_cos):
    return masst_query_all(query_df, database, masst_type, analog, precursor_mz_tol, fragment_mz_tol, min_cos)


def load_and_merge_data(fastmasst_results: pd.DataFrame, usis_table: pd.DataFrame, df_redu: pd.DataFrame,
                        tolerance: float):
    """Load and create the merged dataframe with tolerance filtering"""

    # Early exit if no data
    if len(fastmasst_results) == 0:
        return pd.DataFrame()

    usi_to_name = dict(zip(usis_table['usi'], usis_table['compound_name']))
    compound_names = list(usis_table['compound_name'])

    df_combined = fastmasst_results.copy()
    # filter for tolerance
    df_combined = df_combined[(df_combined['Delta Mass'] >= -tolerance) & (df_combined['Delta Mass'] <= tolerance)]
    df_combined['Compound'] = df_combined['query_usi'].apply(lambda x: usi_to_name.get(x, 'Unknown'))

    # Create filepath column
    usi_parts = df_combined['USI'].str.rsplit('/', n=1).str[-1].str.split(':', n=1).str[0]
    df_combined['filepath'] = df_combined['Dataset'] + "/" + usi_parts.str.replace(
        r'\.mz(ML|XML)$', '', regex=True
    )

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

        df_merged['UBERONBodyPartName'] = df_merged['UBERONBodyPartName'].replace(body_part_replacements)

    return df_merged


@st.cache_data
def load_and_merge_data_wrapper(fastmasst_results: pd.DataFrame, usis_table: pd.DataFrame, df_redu: pd.DataFrame,
                               tolerance: float):
    return load_and_merge_data(fastmasst_results, usis_table, df_redu, tolerance)


# Helper functions
def analyze_counts(df, column_interest):
    """Prepare a table with counts of the fastMASST results."""
    list_body_parts = df[column_interest].unique().tolist()
    df_body_parts = pd.DataFrame(list_body_parts, columns=[column_interest])

    df_counts = df[column_interest].value_counts().rename_axis(column_interest).reset_index(name='Counts_fastMASST')

    compounds = df.groupby(column_interest).agg({
        'Compound': ['nunique', lambda x: ', '.join(map(str, x.unique()))],
        'ATTRIBUTE_DatasetAccession': lambda x: '; '.join(map(str, x.unique()))
    }).reset_index()
    # Flatten MultiIndex columns
    compounds.columns = [column_interest, 'Compounds', 'CompoundsList', 'DatasetAccessions']


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


def filter_by_organism(df_merged: pd.DataFrame, df_redu: pd.DataFrame, organism_filter: str | list = "All organisms"):
    """Filter the merged dataframe according to organism inputs"""

    taxonomy_filters = {
        "Humans": ['9606|Homo sapiens'],
        "Rodents": ['10088|Mus', '10090|Mus musculus', '10105|Mus minutoides', '10114|Rattus',
                    '10116|Rattus norvegicus']
    }

    if isinstance(organism_filter, str):
        if organism_filter in taxonomy_filters:
            selected_taxa = taxonomy_filters[organism_filter]
            df_filtered = df_merged[df_merged['NCBITaxonomy'].isin(selected_taxa)]
            df_redu_filtered = df_redu[df_redu['NCBITaxonomy'].isin(selected_taxa)]
            return df_filtered, df_redu_filtered
        elif organism_filter == "All organisms":
            return df_merged, df_redu
    else:
        # organism_filter is a list
        df_filtered = df_merged[df_merged['NCBITaxonomy'].isin(organism_filter)]
        df_redu_filtered = df_redu[df_redu['NCBITaxonomy'].isin(organism_filter)]
        return df_filtered, df_redu_filtered


def create_heatmap(pivot_table, variable, metric, log_scale=False, normalize=False,
                   col_cluster=False, row_cluster=False, width=2, height=8, custom_sizing=False):
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
        cbar_kws={'orientation': 'vertical'},
        figsize=(width, height) if custom_sizing else None
    )

    fig.ax_heatmap.yaxis.set_ticks_position('left')
    fig.ax_heatmap.yaxis.set_label_position('left')
    fig.ax_heatmap.set_xticklabels(fig.ax_heatmap.get_xticklabels(), rotation=90)
    fig.ax_heatmap.set_yticklabels(fig.ax_heatmap.get_yticklabels(), rotation=0)

    heatmap_title = f"{variable}"

    if col_cluster:
        dendrogram_height = fig.ax_col_dendrogram.get_position().height
        fig.ax_heatmap.set_title(heatmap_title, pad=dendrogram_height * fig.figure.get_dpi() * 5, fontsize=14, weight='bold')
    else:
        fig.ax_heatmap.set_title(heatmap_title, pad=10, fontsize=14, weight='bold')

    if row_cluster:
        dendrogram_width = fig.ax_row_dendrogram.get_position().width
        pad_value = int(dendrogram_width * fig.figure.get_dpi() * 4)  # Convert to points
        print(pad_value)
    else:  
        pad_value = 10  

    fig.ax_heatmap.tick_params(axis='y', pad=pad_value)

    cax = fig.ax_cbar
    # Position colorbar to align with right edge of heatmap, avoiding overlap
    heatmap_pos = fig.ax_heatmap.get_position()
    cax.set_position([heatmap_pos.x1 + 0.04, heatmap_pos.y0 + 0.2, 0.06, heatmap_pos.height * 0.5])
    cbar = fig.ax_heatmap.collections[0].colorbar

    if log_scale:
        cbar.set_label('log2(spectral matches)', fontsize=10)
    elif normalize:
        cbar.set_label('Normalized spectral matches (%)', fontsize=10)
    else:
        cbar.set_label('Spectral matches', fontsize=10)

    fig.ax_heatmap.set(xlabel=None)

    line_count = len(pivot_table[variable].unique())
    fig.ax_heatmap.axhline(y=line_count, color='black', linewidth=1.5)
    fig.ax_heatmap.axvline(x=0, color='black', linewidth=1.5)

    return fig

@st.fragment
def create_heatmap_section(df_filtered, df_redu_filtered, analysis_column):
    """Generate heatmap visualization section with options and controls"""
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
        "Correlation": "correlation",
        "Cosine": "cosine",
        "Euclidean": "euclidean",
        "Jaccard": "jaccard"
    }

    with col2:
        metric = st.selectbox(
            "Clustering metric:",
            options=list(clustering_metrics.keys()),
            help='The distance metric to use. See scipy.spatial.distance.pdist documentation for all available metrics.',
            accept_new_options=True,
        )

    col1, col2 = st.columns(2)
    with col1:
        col_cluster = st.checkbox("Cluster columns", value=False)
    with col2:
        row_cluster = st.checkbox("Cluster rows", value=False)

    col_checkbox, col1, col2 = st.columns([1,1,1])
    with col_checkbox:
        # add space
        st.write("Custom heatmap size")
        custommize_size = st.toggle("Enable", value=False, key='heatmap_size_customize'),
    with col1:
        heatmap_width = st.number_input("Heatmap width", min_value=2, max_value=10, value=4, disabled=not custommize_size)
    with col2:
        heatmap_height = st.number_input("Heatmap height", min_value=4, max_value=15, value=4, disabled=not custommize_size )

    if len(df_filtered) > 0:
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

                if metric in clustering_metrics.keys():
                    metric = clustering_metrics[metric]
                elif metric in ssd._METRICS.keys():
                    metric = metric
                else:
                    metric = 'euclidean'
                    st.toast('Clustering metric not recognized, using default: euclidean', icon="⚠️")

                fig = create_heatmap(
                    pivot_table,
                    analysis_column,
                    metric=metric,
                    log_scale=log_transform,
                    normalize=normalize_redu,
                    col_cluster=col_cluster,
                    row_cluster=row_cluster,
                    width=heatmap_width,
                    height=heatmap_height,
                    custom_sizing=st.session_state.get('heatmap_size_customize', False)
                )
                with st.container(border=True):
                    _, fig_col, _ = st.columns([1, 3, 1])
                    with fig_col:
                        st.pyplot(fig.figure)

                # Download options
                col_a, col_b = st.columns(2)

                with col_a:
                    # Download pivot table
                    csv_pivot = pivot_table.to_csv(index=False)
                    st.download_button(
                        label="Download pivot table",
                        data=csv_pivot,
                        file_name=f"pivot_table_{analysis_column}_{heatmap_type.replace(' ', '_').lower()}.csv",
                        mime="text/csv",
                        icon=":material/download:"
                    )

                with col_b:
                    # Download plot as PNG
                    img_buffer = io.BytesIO()
                    fig.figure.savefig(img_buffer, format='svg', bbox_inches='tight')
                    img_buffer.seek(0)

                    st.download_button(
                        label="Download heatmap (SVG)",
                        data=img_buffer.getvalue(),
                        file_name=f"heatmap_{analysis_column}_{heatmap_type.replace(' ', '_').lower()}.svg",
                        mime="image/svg+xml",
                        icon=":material/download:",
                    )

            except Exception as e:
                st.error(f"Error generating heatmap: {str(e)}")
    else:
        st.info("No data available to generate heatmap with current filters.", icon="ℹ️")


@st.fragment
def create_bar_chart_section(df_filtered, df_redu_filtered, analysis_column, compound_name):
    """Generate bar chart visualization for single USI queries"""
    import matplotlib.pyplot as plt
    
    st.header("📊 Bar Chart Visualization")
    st.caption(f"Distribution of spectral matches for **{compound_name}**")

    # Bar chart options
    col1, col2 = st.columns(2)
    with col1:
        chart_type = st.selectbox(
            "Chart type:",
            options=['Raw counts', 'ReDU-normalized counts'],
            help="Choose the type of bar chart to generate",
            key='bar_chart_type'
        )
    
    with col2:
        sort_order = st.selectbox(
            "Sort by:",
            options=['Descending (highest first)', 'Ascending (lowest first)', 'Alphabetical'],
            help="Choose how to sort the bars",
            key='bar_sort_order'
        )

    col1, col2 = st.columns(2)
    with col1:
        orientation = st.radio(
            "Orientation:",
            options=['Horizontal', 'Vertical'],
            horizontal=True,
            key='bar_orientation'
        )
    with col2:
        show_values = st.checkbox("Show values on bars", value=True, key='bar_show_values')

    if len(df_filtered) > 0:
        with st.spinner("Generating bar chart..."):
            try:
                # Apply exclusions if any
                if variable_to_exclude:
                    drop_mask = df_filtered[analysis_column].isin(variable_to_exclude)
                    indices_to_drop = df_filtered[drop_mask].index
                    df_filtered = df_filtered.drop(indices_to_drop)

                # Get counts
                counts_data = df_filtered[analysis_column].value_counts().reset_index()
                counts_data.columns = [analysis_column, 'Counts']

                # Apply ReDU normalization if selected
                normalize_redu = chart_type == 'ReDU-normalized counts'
                if normalize_redu:
                    df_redu_counts = df_redu_filtered[analysis_column].value_counts().reset_index()
                    df_redu_counts.columns = [analysis_column, 'ReDU_Counts']
                    
                    counts_data = pd.merge(counts_data, df_redu_counts, on=analysis_column, how='left')
                    counts_data['Normalized'] = (counts_data['Counts'] / counts_data['ReDU_Counts']) * 100
                    counts_data = counts_data.dropna(subset=['Normalized'])
                    value_column = 'Normalized'
                    ylabel = 'Normalized spectral matches (%)'
                else:
                    value_column = 'Counts'
                    ylabel = 'Spectral matches'

                # Sort data
                if sort_order == 'Descending (highest first)':
                    counts_data = counts_data.sort_values(value_column, ascending=False)
                elif sort_order == 'Ascending (lowest first)':
                    counts_data = counts_data.sort_values(value_column, ascending=True)
                else:  # Alphabetical
                    counts_data = counts_data.sort_values(analysis_column, ascending=True)

                # Create figure
                num_categories = len(counts_data)
                if orientation == 'Horizontal':
                    fig_height = max(4, num_categories * 0.4)
                    fig_width = 10
                else:
                    fig_width = max(8, num_categories * 0.5)
                    fig_height = 6

                fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                # Create bar plot
                if orientation == 'Horizontal':
                    bars = ax.barh(
                        counts_data[analysis_column],
                        counts_data[value_column],
                        color=sns.color_palette("Blues_r", n_colors=num_categories),
                        edgecolor='white',
                        linewidth=0.5
                    )
                    ax.set_xlabel(ylabel)
                    ax.set_ylabel(analysis_column)
                    ax.invert_yaxis()  # Highest at top
                    
                    if show_values:
                        for bar, val in zip(bars, counts_data[value_column]):
                            width = bar.get_width()
                            label = f'{val:.1f}%' if normalize_redu else f'{int(val)}'
                            ax.text(width + (ax.get_xlim()[1] * 0.01), bar.get_y() + bar.get_height()/2,
                                   label, va='center', ha='left', fontsize=9)
                else:
                    bars = ax.bar(
                        counts_data[analysis_column],
                        counts_data[value_column],
                        color=sns.color_palette("Blues_r", n_colors=num_categories),
                        edgecolor='white',
                        linewidth=0.5
                    )
                    ax.set_ylabel(ylabel)
                    ax.set_xlabel(analysis_column)
                    plt.xticks(rotation=45, ha='right')
                    
                    if show_values:
                        for bar, val in zip(bars, counts_data[value_column]):
                            height = bar.get_height()
                            label = f'{val:.1f}%' if normalize_redu else f'{int(val)}'
                            ax.text(bar.get_x() + bar.get_width()/2, height + (ax.get_ylim()[1] * 0.01),
                                   label, ha='center', va='bottom', fontsize=9, rotation=0)

                ax.set_title(f'{compound_name} - Distribution by {analysis_column}', fontsize=14, weight='bold', pad=10)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                
                plt.tight_layout()

                # Display chart
                with st.container(border=True):
                    st.pyplot(fig)

                # Download options
                col_a, col_b = st.columns(2)

                with col_a:
                    # Download data table
                    csv_data = counts_data.to_csv(index=False)
                    st.download_button(
                        label="Download data table",
                        data=csv_data,
                        file_name=f"bar_chart_data_{analysis_column}_{chart_type.replace(' ', '_').lower()}.csv",
                        mime="text/csv",
                        icon=":material/download:"
                    )

                with col_b:
                    # Download plot as SVG
                    img_buffer = io.BytesIO()
                    fig.savefig(img_buffer, format='svg', bbox_inches='tight')
                    img_buffer.seek(0)

                    st.download_button(
                        label="Download chart (SVG)",
                        data=img_buffer.getvalue(),
                        file_name=f"bar_chart_{analysis_column}_{chart_type.replace(' ', '_').lower()}.svg",
                        mime="image/svg+xml",
                        icon=":material/download:",
                    )

                plt.close(fig)

            except Exception as e:
                st.error(f"Error generating bar chart: {str(e)}")
    else:
        st.info("No data available to generate bar chart with current filters.", icon="ℹ️")


@st.fragment
def render_datasets_distribution(df_filtered, analisys_col, unique_organs):
    """Render dataset distribution bar chart"""

    if len(df_filtered) > 0:
        col1, col2 = st.columns(2)
        with col1:
            compounds = df_filtered['Compound'].unique().tolist()
            compounds.insert(0, "All compounds")
            compound_selection = st.selectbox("Select compound", compounds, key="dataset_distribution_compound")
        with col2:
            unique_options = df_filtered[analisys_col].unique().tolist()
            unique_options.sort()
            unique_options.insert(0, "All")
            strata_selection = st.selectbox("Stratify by", unique_options, key="stratify_by")

        with st.spinner("Generating dataset distribution chart..."):
            try:
                # Apply filters based on selections
                dataset_subset = df_filtered.copy()
                
                if compound_selection != "All compounds":
                    dataset_subset = dataset_subset[dataset_subset['Compound'] == compound_selection]
                
                if strata_selection != "All":
                    dataset_subset = dataset_subset[dataset_subset[analisys_col] == strata_selection]
                
                dataset_counts = dataset_subset['Dataset'].value_counts().reset_index()
                dataset_counts.columns = ['Dataset', 'Counts']

                fig, ax = plt.subplots(figsize=(10, 6))
                bars = ax.barh(
                    dataset_counts['Dataset'],
                    dataset_counts['Counts'],
                    color=sns.color_palette("Greens_r", n_colors=len(dataset_counts)),
                    edgecolor='white',
                    linewidth=0.5
                )
                ax.set_xlabel(f'Spectral matches for {compound_selection}')
                ax.set_ylabel('Dataset')
                ax.invert_yaxis()  # Highest at top

                for bar, val in zip(bars, dataset_counts['Counts']):
                    width = bar.get_width()
                    ax.text(width + (ax.get_xlim()[1] * 0.01), bar.get_y() + bar.get_height()/2,
                            f'{int(val)}', va='center', ha='left', fontsize=9)

                # Build title based on selections
                title_parts = []
                if compound_selection != "All compounds":
                    title_parts.append(compound_selection)
                if strata_selection != "All":
                    title_parts.append(strata_selection)
                chart_title = ' - '.join(title_parts) if title_parts else 'All Data'
                
                ax.set_title(f'Spectral Matches: {chart_title}', fontsize=14, weight='bold', pad=10)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)

                plt.tight_layout()

                with st.container(border=True):
                    st.pyplot(fig)

                # Download options
                csv_data = dataset_counts.to_csv(index=False)
                st.download_button(
                    label="Download dataset distribution data",
                    data=csv_data,
                    file_name="dataset_distribution_data.csv",
                    mime="text/csv",
                    icon=":material/download:"
                )

                img_buffer = io.BytesIO()
                fig.savefig(img_buffer, format='svg', bbox_inches='tight')
                img_buffer.seek(0)

                st.download_button(
                    label="Download dataset distribution chart (SVG)",
                    data=img_buffer.getvalue(),
                    file_name="dataset_distribution_chart.svg",
                    mime="image/svg+xml",
                    icon=":material/download:",
                )

                plt.close(fig)

            except Exception as e:
                st.error(f"Error generating dataset distribution chart: {str(e)}")
    else:
        st.info("No data available to generate dataset distribution chart with current filters.", icon="ℹ️")


@st.fragment
def render_dataset_details_card(unique_datasets):
    """Render dataset details card"""
    st.header("🗂️ Dataset Details")
    dataset_id = st.selectbox("Select datasets to view details:", options=unique_datasets, key="massive_dataset_selection")

    
    if dataset_id:
        display_metabolomics_metadata(dataset_id)
    else:

        st.markdown(
            """
            The datasets analyzed in this tool are sourced from the [MASSIVE](https://massive.ucsd.edu/) repository, which is a comprehensive resource for mass spectrometry data. 
            For more information about a specific dataset, please visit the MASSIVE website and use the dataset accession number provided in the results.
            
            **How to find dataset details:**
            1. Go to the [MASSIVE website](https://massive.ucsd.edu/).
            2. Use the search bar to enter the dataset accession number (e.g., MSV000083437).
            3. Explore the dataset page for detailed information including experimental design, sample metadata, and associated publications.
            
            For any questions or further assistance, please refer to the [MASSIVE FAQ](https://massive.ucsd.edu/FAQ.jsp) or contact the MASSIVE support team.
            """
        )

# Sidebar for inputs
# Load query parameters from URL (only on first load)
if '_query_params_loaded' not in st.session_state:
    st.session_state._url_params = load_all_params()
    st.session_state._query_params_loaded = True
else:
    st.session_state._url_params = st.session_state.get('_url_params', {})

url_params = st.session_state._url_params

with st.sidebar:
    st.header("📁 Data Input")

    # Determine initial example state from URL params
    initial_example = url_params.get('example')
    
    # Example data option
    use_example = st.checkbox("Load example data", 
                              value=initial_example is not None, 
                              help="Use built-in example files instead of uploading your own", 
                              key='use_example'
                              )

    if use_example:

        def clean_example_data():
            st.session_state.pop("results", None)
            st.session_state.pop("df_merged", None)

        # Set initial selection based on URL param
        example_options = ["Example Set 1", "Example Set 2"]
        if initial_example == '2':
            initial_example_idx = 1
        else:
            initial_example_idx = 0
            
        selected_example = st.selectbox(
            "Select example dataset:",
            options=example_options,
            index=initial_example_idx,
            help="Choose which example dataset to load",
            on_change=clean_example_data,
        )
        if selected_example == "Example Set 1":
            data = pd.DataFrame(EXAMPLE_SET_1)
            # Sync example selection to URL
            sync_param('example', '1')
        elif selected_example == 'Example Set 2':
            data = pd.DataFrame(EXAMPLE_SET_2)
            sync_param('example', '2')
    else:
        # Load USI data from URL params if available
        data = url_params.get('usi')
        # Clear example from URL when not using example
        if 'example' in st.query_params:
            del st.query_params['example']

    usi_data = create_usi_input(usi_data=data, sync_url=not use_example)
    masst_query_params = create_masst_sidebar(initial_params=url_params)
    # Filter out empty rows
    query_df = usi_data[usi_data['usi'].str.strip() != ''].copy()

    st.session_state.masst_query_params = masst_query_params

    # Change button label and behavior based on example mode
    if use_example:
        button_label = "Load Example Data"
        button_icon = "🧪"
    else:
        button_label = "Run MASST Query"
        button_icon = "🔎"

    st.session_state['run_masst_query'] = st.button(
        button_label, 
        icon=button_icon, 
        type="primary", 
        use_container_width=True,
        disabled='results' in st.session_state  # Disable if example already loaded
    )

    if st.button("Run new query", type="secondary", icon=":material/replay:",use_container_width=True):
        st.session_state.clear()
        st.session_state["use_example"] = False
        # Clear URL query parameters
        clear_all_params()
        st.rerun()

    st.markdown(f"**ReDU metadata last updated:** {file_date}")

    st.subheader("Contributors")
    st.markdown(
        """
    - [Vincent Lamoureux PhD](https://scholar.google.com/citations?user=_OboZ0YAAAAJ) - UC San Diego
    - [Helena Russo PhD](https://sites.google.com/view/helenamrusso/home) - UC San Diego
    - [Prajit Rajkumar](https://scholar.google.com/citations?user=_iKPhb0AAAAJ) - UC San Diego
    - [Wilhan Nunes PhD](https://scholar.google.com/citations?user=4cPVoeIAAAAJ) - UC San Diego
    - [Mingxun Wang PhD](https://www.cs.ucr.edu/~mingxunw/) - UC Riverside
    """
    )

    st.subheader("Documentations and Resources")
    st.markdown("""
        - [Other Similar Tools](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/toolindex/)
        - [MASST Documentation](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/masst/)
        - [Fast Search Page](https://fasst.gnps2.org/fastsearch/)
        """)

    with st.expander("API Reference"):
        st.markdown("**Parameter Schema**")
        st.json(get_param_schema())


if st.session_state.get('run_masst_query', False):
    # Skip actual MASST query if using example data (already loaded)
    if use_example:
        results = cached_masst_results(query_df, **masst_query_params)
        st.session_state.results = results

    elif 0 < len(query_df) <= 20:
        with st.spinner("Running MASST query..."):
            # Here you would call your imported masst_query_all function
            results = masst_query_all(query_df, **masst_query_params)
            st.success(f"Query performed with {len(query_df)} USIs with parameters: {masst_query_params}")
            st.session_state.results = results
            # st.write(results)
    elif len(query_df) > 20:
        st.error("The query is limited to 20 USIs. Please, use the [Fast Batch](https://gnps2.org/workflowinput?workflowname=fasst_batch_workflow)"
                 " workflow fore more flexibility.", icon="⚠️")
    else:
        st.error("Please add at least one USI to query", icon="⚠️")

if "results" in st.session_state and len(st.session_state.results) > 0:
    # Use df_redu from session state
    if df_redu is None:
        st.error("Failed to load ReDU metadata file")
        st.stop()

    results = st.session_state.results

    try:
        st.markdown("----")
        # Organism selection
        col1, col2 = st.columns(2)
        with col1:
            st.write("🧬 Organism Filter")

            # available_organisms = df_merged['NCBITaxonomy'].dropna().unique()

            # Get initial organism from URL params
            initial_organism = url_params.get('organism', get_default('organism'))
            organism_options = ['Humans', 'Rodents', 'All organisms', "Manual Entry"]
            if initial_organism in organism_options:
                initial_org_idx = organism_options.index(initial_organism)
            else:
                initial_org_idx = 0  # Default to Humans

            organism_choice = st.selectbox(
                "Select organism:",
                options=organism_options,
                index=initial_org_idx,
                help="Filter data by organism type",
                key='organism_choice'
            )
            
            # Sync organism to URL
            sync_param('organism', organism_choice)

            if organism_choice == "Manual Entry":
                organism_choice = st.multiselect(
                    "Search NCBI Taxonomy ID or scientific to filter the list",
                    options=df_redu['NCBITaxonomy'].dropna().unique().tolist(),
                    default=['9606|Homo sapiens'],
                    help="Enter a specific NCBI Taxonomy ID or scientific name (multiselect)",
                    key='manual_organism_choice',
                )
        if "df_merged" not in st.session_state:
            start_time = time.time()
            
            # Use cached version for example data, uncached for real queries
            if use_example:
                df_merged = load_and_merge_data_wrapper(results, usi_data, df_redu, masst_query_params.get('precursor_mz_tol', 0.02))
            else:
                df_merged = load_and_merge_data(results, usi_data, df_redu, masst_query_params.get('precursor_mz_tol', 0.02))
                end_time = time.time()
                elapsed_time = end_time - start_time
                st.success(f"✅ Data merged in {elapsed_time:.2f} seconds")
            
            # save for example
            # df_merged.to_csv('example_data/df_merged_example.csv', index=False)
            st.session_state.df_merged = df_merged
        else:
            df_merged = st.session_state.df_merged

        df_filtered, df_redu_filtered = filter_by_organism(df_merged, df_redu, organism_filter=organism_choice)
        rename_dict = {"UBERONBodyPartName": "Body Part",
                       "DOIDCommonName": "Disease"
                       }
        df_filtered.rename(columns=rename_dict, inplace=True)
        df_redu_filtered.rename(columns=rename_dict, inplace=True)

        with col2:
            # Analysis options
            st.write("📈 Analysis Options")

            available_columns = [col for col in df_filtered.columns if
                                 col in ['Body Part', 'Disease', 'ATTRIBUTE_SubjectGender',
                                         'ATTRIBUTE_Age']]

            analysis_column = st.selectbox(
                "Select analysis variable:",
                options=available_columns,
                help="Choose which ReDU column to analyze"
            )

            unique_options = list(df_filtered[analysis_column].unique())
            variable_to_exclude = st.multiselect("Exclude from analysis (optional)", options=unique_options)
        if len(df_filtered) == 0:
            st.warning(f"No data available for the selected organism filter: {', '.join(organism_choice)}. "
                       f"Try selecting 'All organisms' or a different filter.", icon="⚠️")
        else:
            if isinstance(organism_choice, list):
                st.info(f"📊 Filtered data: {len(df_filtered):,} spectral matches for {', '.join(organism_choice)}")
            else:
                st.info(f"📊 Filtered data: {len(df_filtered):,} spectral matches for {organism_choice}")

        # Main content
        col1, col2 = st.columns([1, 2])

        with col1:
            st.header("📋 Summary")

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
                    mime="text/csv",
                    disabled= not len(counts_table) > 0,
                )


        with col2:
            # Determine if single or multiple USI query
            num_compounds = df_filtered['Compound'].nunique()
            
            if num_compounds == 1:
                # Single USI - show bar chart
                compound_name = df_filtered['Compound'].iloc[0]
                create_bar_chart_section(df_filtered, df_redu_filtered, analysis_column, compound_name)
            else:
                # Multiple USIs - show heatmap
                create_heatmap_section(df_filtered, df_redu_filtered, analysis_column)

        # Data preview
        st.header(":mag: Data Preview")
        preview_col, datasets_fig_col = st.columns([1, 1])
        with datasets_fig_col:
            # Datasets distribution plot
            st.subheader("Datasets Distribution", help="Distribution of spectral matches across datasets")
            render_datasets_distribution(df_filtered, analysis_column, unique_options)
        with preview_col:
            st.subheader("Merged Dataset", help="Showing first 100 rows of the merged dataset")
            st.dataframe(df_filtered.head(100))

            # Download full dataset
            csv_full = df_filtered.to_csv(index=False)
            st.download_button(
                label="📥 Download full merged dataset",
                data=csv_full,
                file_name=f"merged_dataset.csv",
                mime="text/csv",
                disabled= not len(df_filtered) > 0,
            )
        
        unique_datasets = df_filtered['Dataset'].unique().tolist()
        render_dataset_details_card(unique_datasets)

    except Exception as e:
        st.error(f"❌ Error processing data: {str(e)}")
        st.info("Please check your file formats and ensure they match the expected structure.")
        raise

elif "results" in st.session_state and len(st.session_state.get('results', [])) == 0:
    st.warning("MASST query returned no results. Try again modifying the parameters in the sidebar or using a different set of USIs.", icon="⚠️")

else:
    st.info("👈 Please edit the USI input data table in the sidebar, choose your parameters and click Run MASST Query to start.")

    welcome_page()
