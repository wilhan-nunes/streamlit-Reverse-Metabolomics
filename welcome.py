import streamlit as st


def welcome_page():
    st.markdown("""
    ### 📖 Citation
    Charron-Lamoureux, V., Mannochio-Russo, H., Lamichhane, S. et al. A guide to reverse metabolomics—a framework for big data discovery strategy. Nat Protoc (2025). https://doi.org/10.1038/s41596-024-01136-2
    
    ### 📘 How to use this tool:

    1. **Enter or upload USIs**: In the sidebar, input your USIs (Universal Spectrum Identifiers) or load the example data.
    2. **Set MASST parameters**: Adjust search parameters as needed.
    3. **Run MASST Query**: Click the "Run MASST Query" button to retrieve fastMASST results.
    4. **Adjust analysis settings**: Set delta mass tolerance, select organism, and choose the analysis variable (e.g., body part, disease, gender, or age).
    5. **Generate visualizations**: Select heatmap options and click "Generate Heatmap" to view interactive results.
    6. **Download results**: Export summary tables and heatmaps as CSV or PNG for further analysis.

    ### 📊 Available visualizations:
    - **Raw counts**: Direct spectral match counts.
    - **Log-transformed counts**: Log2-transformed counts for improved visualization of low-abundance features.
    - **ReDU-normalized counts**: Counts normalized by sample availability in the ReDU database.

    ### 💡 Tips:
    - Use clustering options to reveal patterns in your data.
    - Download tables and images for offline analysis or publication.
    - Example data is available for quick testing and demonstration.
    """)
