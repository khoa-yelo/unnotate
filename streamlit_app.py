"""
Streamlit version of the protein accession visualizer
"""
import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import os
import zipfile
import tempfile
import io
from pathlib import Path
from matplotlib import cm, colors as mcolors

# Helper for categorical color mapping
def get_categorical_colormap(categories):
    base_colors = list(mcolors.TABLEAU_COLORS.values()) + list(mcolors.CSS4_COLORS.values())
    color_map = {}
    for i, cat in enumerate(categories):
        color_map[cat] = base_colors[i % len(base_colors)]
    return color_map

def process_uploaded_zip(uploaded_file):
    """Process uploaded zip file and extract data"""
    try:
        with zipfile.ZipFile(uploaded_file, 'r') as zip_ref:
            # List all files in the zip
            file_list = zip_ref.namelist()
            st.info(f"Files found in zip: {file_list}")
            
            # Extract files to temporary directory
            with tempfile.TemporaryDirectory() as temp_dir:
                zip_ref.extractall(temp_dir)
                
                # Look for the required files
                csv_file = None
                accession_file = None
                similarity_file = None
                
                for file_path in file_list:
                    if file_path.endswith('.csv'):
                        csv_file = os.path.join(temp_dir, file_path)
                    elif 'accession' in file_path.lower() and file_path.endswith('.npy'):
                        accession_file = os.path.join(temp_dir, file_path)
                    elif 'similarity' in file_path.lower() and file_path.endswith('.npy'):
                        similarity_file = os.path.join(temp_dir, file_path)
                
                if not csv_file:
                    st.error("❌ No CSV file found in the zip")
                    return None, None, None
                
                if not accession_file:
                    st.error("❌ No accession arrays file found in the zip")
                    return None, None, None
                
                if not similarity_file:
                    st.error("❌ No similarity array file found in the zip")
                    return None, None, None
                
                # Load the data
                df = pd.read_csv(csv_file)
                st.success(f"✅ Loaded CSV with {len(df)} rows")
                
                accession_arrays_loaded = np.load(accession_file)
                similarity_array = np.load(similarity_file)
                
                # Convert to list format
                accession_arrays = []
                for arr in accession_arrays_loaded:
                    if hasattr(arr, 'tolist'):
                        arr = arr.tolist()
                    arr = [acc.decode('utf-8') if isinstance(acc, bytes) else acc for acc in arr]
                    accession_arrays.append(arr)
                
                st.success(f"✅ Loaded {len(accession_arrays)} accession arrays and similarity array")
                
                return df, accession_arrays, similarity_array
                
    except Exception as e:
        st.error(f"❌ Error processing zip file: {e}")
        return None, None, None

@st.cache_data
def load_data(dataset_type="Virus"):
    """Load and cache the protein data"""
    # Set up path logic for cloud deployment
    BASE_DIR = Path(__file__).resolve().parent  
    data_path = BASE_DIR / "data"
    
    # Determine file prefix based on dataset type
    prefix = "viral" if dataset_type == "Virus" else "bacterial"
    
    # Try multiple possible paths for the CSV file
    csv_paths = [
        f"{prefix}_parsed_uniprot_swiss_data.csv",
        f"data/{prefix}_parsed_uniprot_swiss_data.csv",
        str(data_path / f"{prefix}_parsed_uniprot_swiss_data.csv"),
        f"../data/{prefix}_parsed_uniprot_swiss_data.csv",
        f"data/unnotate/{prefix}_parsed_uniprot_swiss_data.csv",
        f"../data/unnotate/{prefix}_parsed_uniprot_swiss_data.csv",
    ]
    
    df = None
    for path in csv_paths:
        if os.path.exists(path):
            df = pd.read_csv(path)
            break
    
    if df is None:
        return None, None, None
    
    # Try to load numpy arrays
    accession_arrays = None
    similarity_array = None
    
    # Try multiple possible paths for numpy files
    npy_paths = [
        (f"{prefix}_accession_arrays.npy", f"{prefix}_similarity_array.npy"),
        (f"data/{prefix}_accession_arrays.npy", f"data/{prefix}_similarity_array.npy"),
        (str(data_path / f"{prefix}_accession_arrays.npy"), str(data_path / f"{prefix}_similarity_array.npy")),
        (f"../data/{prefix}_accession_arrays.npy", f"../data/{prefix}_similarity_array.npy"),
    ]
    
    for acc_path, sim_path in npy_paths:
        if os.path.exists(acc_path) and os.path.exists(sim_path):
            try:
                accession_arrays_loaded = np.load(acc_path)
                similarity_array = np.load(sim_path)
                
                # Convert to list format
                accession_arrays = []
                for arr in accession_arrays_loaded:
                    if hasattr(arr, 'tolist'):
                        arr = arr.tolist()
                    arr = [acc.decode('utf-8') if isinstance(acc, bytes) else acc for acc in arr]
                    accession_arrays.append(arr)
                
                break
                
            except Exception as e:
                continue
    
    if accession_arrays is None:
        return None, None, None
    
    return df, accession_arrays, similarity_array

def create_heatmap(df, accession_arrays, similarity_array, selected_domain=None):
    """Create heatmap visualization with optional domain-based sorting"""
    # Flatten all accessions for batch processing first
    all_accessions = []
    for accessions in accession_arrays:
        all_accessions.extend(accessions)
    
    # Batch lookup using pandas isin() - do this once
    mask = df['accession'].isin(all_accessions)
    relevant_df = df[mask].set_index('accession')
    
    # If a domain is selected, sort the arrays
    if selected_domain and selected_domain != "All Domains":
        sorted_arrays = []
        sorted_similarities = []
        
        for i, accessions in enumerate(accession_arrays):
            similarities = similarity_array[i].tolist() if hasattr(similarity_array[i], 'tolist') else similarity_array[i]
            
            # Get protein data for sorting using vectorized lookup
            protein_data = []
            for j, acc in enumerate(accessions):
                if acc in relevant_df.index:
                    row = relevant_df.loc[acc]
                    if pd.notna(row['taxonomy_lineage']):
                        taxonomy = row['taxonomy_lineage']
                        domain = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
                    else:
                        domain = 'Unknown'
                else:
                    domain = 'Unknown'
                protein_data.append((acc, domain, similarities[j]))
            
            # Sort by domain first, then by similarity
            def sort_key(item):
                acc, domain, sim = item
                if domain == selected_domain:
                    return (0, -sim)  # Selected domain first, then by similarity (descending)
                else:
                    return (1, domain, -sim)  # Other domains, sorted by domain name, then similarity
            
            sorted_data = sorted(protein_data, key=sort_key)
            sorted_accessions = [acc for acc, domain, sim in sorted_data]
            sorted_sims = [sim for acc, domain, sim in sorted_data]
            
            sorted_arrays.append(sorted_accessions)
            sorted_similarities.append(sorted_sims)
        
        # Use sorted arrays for visualization
        arrays_to_visualize = sorted_arrays
    else:
        # Use original arrays
        arrays_to_visualize = accession_arrays
    
    max_len = max(len(arr) for arr in arrays_to_visualize)
    
    # Create domain data efficiently
    domain_data = []
    hover_texts = []
    customdata = []
    
    for cds_idx, accessions in enumerate(arrays_to_visualize):
        row_domains = []
        row_hover = []
        row_urls = []
        
        for pos_idx, acc in enumerate(accessions):
            # Use vectorized lookup
            if acc in relevant_df.index:
                row = relevant_df.loc[acc]
                if pd.notna(row['taxonomy_lineage']):
                    taxonomy = row['taxonomy_lineage']
                    domain = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
                else:
                    domain = 'Unknown'
                
                # Create hover text efficiently
                function_text = row['function']
                if pd.notna(function_text):
                    function_display = function_text[:50] + ('...' if len(function_text) > 50 else '')
                else:
                    function_display = 'No function data'
                
                hover_text = f"""
                <b>CDS Index:</b> {cds_idx+1}<br>
                <b>Position Index:</b> {pos_idx+1}<br>
                <b>Accession:</b> {acc}<br>
                <b>Name:</b> {row['name']}<br>
                <b>Full Name:</b> {row['full_name'] if pd.notna(row['full_name']) else 'Not available'}<br>
                <b>Organism:</b> {row['organism_common']}<br>
                <b>Domain:</b> {domain}<br>
                <b>Length:</b> {row['sequence_length'] if pd.notna(row['sequence_length']) else 'Not available'}<br>
                <b>Function:</b> {function_display}<br>
                <b>Click to open UniProt page</b>
                """
            else:
                domain = 'Unknown'
                hover_text = f"<b>Accession:</b> {acc}<br><b>Status:</b> Not found in database<br><b>Click to open UniProt page</b>"
            
            row_domains.append(domain)
            row_urls.append(f"https://www.uniprot.org/uniprot/{acc}")
            row_hover.append(hover_text)
        
        # Pad with empty strings
        row_domains.extend([''] * (max_len - len(accessions)))
        row_hover.extend([''] * (max_len - len(accessions)))
        row_urls.extend([''] * (max_len - len(accessions)))
        domain_data.append(row_domains)
        hover_texts.append(row_hover)
        customdata.append(row_urls)
    
    # Transpose and flip
    domain_data_transposed = list(zip(*domain_data))[::-1]
    hover_texts_transposed = list(zip(*hover_texts))[::-1]
    customdata_transposed = list(zip(*customdata))[::-1]
    
    # Get unique domains efficiently
    unique_domains = sorted(set(d for row in domain_data_transposed for d in row if d and d != ''))
    
    # Create numerical mapping
    domain_to_num = {domain: i for i, domain in enumerate(unique_domains)}
    domain_numerical = []
    for row in domain_data_transposed:
        num_row = []
        for domain in row:
            if domain and domain != '':
                num_row.append(domain_to_num[domain])
            else:
                num_row.append(-1)
        domain_numerical.append(num_row)
    
    # Create color map
    domain_color_map = get_categorical_colormap(unique_domains)
    colorscale = []
    for i, d in enumerate(unique_domains):
        color = domain_color_map[d]
        if len(unique_domains) == 1:
            colorscale = [[0, color], [1, color]]
        else:
            val = i / (len(unique_domains) - 1)
            colorscale.append([val, color])
    if len(unique_domains) > 1:
        last_color = domain_color_map[unique_domains[-1]]
        colorscale.append([1.0, last_color])
    
    # Create figure
    fig = go.Figure(data=go.Heatmap(
        z=domain_numerical,
        hovertext=hover_texts_transposed,
        hoverinfo='text',
        colorscale=colorscale,
        showscale=False,
        zmin=0,
        zmax=len(unique_domains)-1,
        customdata=customdata_transposed
    ))
    
    # Add legend
    for d in unique_domains:
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(size=15, color=domain_color_map[d]),
            legendgroup=d,
            showlegend=True,
            name=d
        ))
    
    title = "Domain Heatmap"
    if selected_domain and selected_domain != "All Domains":
        title += f" - Sorted by {selected_domain} first, then by similarity"
    
    fig.update_layout(
        title=title,
        xaxis_title="CDS Index",
        yaxis_title="Closest match",
        height=400,
        width=1200,
        legend_title="Domain"
    )
    return fig

def create_similarity_heatmap(df, accession_arrays, similarity_array, selected_domain=None):
    """Create similarity heatmap with optional domain-based sorting"""
    # Flatten all accessions for batch processing first
    all_accessions = []
    for accessions in accession_arrays:
        all_accessions.extend(accessions)
    
    # Batch lookup using pandas isin() - do this once
    mask = df['accession'].isin(all_accessions)
    relevant_df = df[mask].set_index('accession')
    
    # If a domain is selected, sort the arrays
    if selected_domain and selected_domain != "All Domains":
        sorted_arrays = []
        sorted_similarities = []
        
        for i, accessions in enumerate(accession_arrays):
            similarities = similarity_array[i].tolist() if hasattr(similarity_array[i], 'tolist') else similarity_array[i]
            
            # Get protein data for sorting using vectorized lookup
            protein_data = []
            for j, acc in enumerate(accessions):
                if acc in relevant_df.index:
                    row = relevant_df.loc[acc]
                    if pd.notna(row['taxonomy_lineage']):
                        taxonomy = row['taxonomy_lineage']
                        domain = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
                    else:
                        domain = 'Unknown'
                else:
                    domain = 'Unknown'
                protein_data.append((acc, domain, similarities[j]))
            
            # Sort by domain first, then by similarity
            def sort_key(item):
                acc, domain, sim = item
                if domain == selected_domain:
                    return (0, -sim)  # Selected domain first, then by similarity (descending)
                else:
                    return (1, domain, -sim)  # Other domains, sorted by domain name, then similarity
            
            sorted_data = sorted(protein_data, key=sort_key)
            sorted_accessions = [acc for acc, domain, sim in sorted_data]
            sorted_sims = [sim for acc, domain, sim in sorted_data]
            
            sorted_arrays.append(sorted_accessions)
            sorted_similarities.append(sorted_sims)
        
        # Use sorted similarities for visualization
        similarity_data = np.array(sorted_similarities)
        arrays_to_visualize = sorted_arrays
    else:
        # Use original similarity array
        similarity_data = similarity_array
        arrays_to_visualize = accession_arrays
    
    # Create customdata with UniProt URLs
    max_len = max(len(arr) for arr in arrays_to_visualize)
    customdata = []
    
    for accessions in arrays_to_visualize:
        row_urls = []
        for acc in accessions:
            uniprot_url = f"https://www.uniprot.org/uniprot/{acc}"
            row_urls.append(uniprot_url)
        
        # Pad with empty strings
        row_urls.extend([''] * (max_len - len(accessions)))
        customdata.append(row_urls)
    
    # Transpose and flip
    customdata_transposed = list(zip(*customdata))[::-1]
    
    z = np.array(similarity_data).T[::-1]
    fig = go.Figure(data=go.Heatmap(
        z=z,
        colorscale='Blues',
        colorbar=dict(title="Similarity"),
        customdata=customdata_transposed
    ))
    
    title = "Similarity Heatmap"
    if selected_domain and selected_domain != "All Domains":
        title += f" - Sorted by {selected_domain} first, then by similarity"
    
    fig.update_layout(
        title=title,
        xaxis_title="Array Index",
        yaxis_title="Position in Array",
        height=400,
        width=1200
    )
    return fig

def main():
    st.set_page_config(page_title="Protein Accession Visualizer", layout="wide")
    
    st.title("Protein Accession Visualization Dashboard")
    
    # Dataset selection dropdown
    st.sidebar.header("Dataset Selection")
    dataset_type = st.sidebar.selectbox(
        "Choose Dataset:",
        ["Virus", "Bacteria"],
        help="Select which dataset to load from local files"
    )
    
    # Always show upload interface
    with st.expander("📥 Upload Data from ZIP File", expanded=True):
        st.write("Upload a zip file containing the required data files:")
        st.write("- A CSV file (e.g., `PREFIX_parsed_uniprot_swiss_data.csv`)")
        st.write("- A numpy file with accession arrays (e.g., `PREFIX_accession_arrays.npy`)")
        st.write("- A numpy file with similarity data (e.g., `PREFIX_similarity_array.npy`)")
        
        uploaded_file = st.file_uploader("Choose a ZIP file", type=["zip"])
        
        if uploaded_file is not None:
            with st.spinner("Processing uploaded ZIP file..."):
                df, accession_arrays, similarity_array = process_uploaded_zip(uploaded_file)
                if df is not None and accession_arrays is not None and similarity_array is not None:
                    st.success("✅ Data loaded successfully from uploaded ZIP file!")
                else:
                    st.error("❌ Failed to load data from uploaded ZIP file. Please check the file contents.")
        else:
            # Try to load from local files if no upload
            # Set up path logic for cloud deployment
            BASE_DIR = Path(__file__).resolve().parent
            data_path = BASE_DIR / "data"
            
            # Determine file prefix based on dataset type
            prefix = "viral" if dataset_type == "Virus" else "bacterial"
            
            # Check multiple possible file locations
            possible_files = [
                # Current directory
                (f"{prefix}_accession_arrays.npy", f"{prefix}_similarity_array.npy", f"{prefix}_parsed_uniprot_swiss_data.csv"),
                # Data directory
                (str(data_path / f"{prefix}_accession_arrays.npy"), str(data_path / f"{prefix}_similarity_array.npy"), str(data_path / f"{prefix}_parsed_uniprot_swiss_data.csv")),
                # Relative paths
                (f"data/{prefix}_accession_arrays.npy", f"data/{prefix}_similarity_array.npy", f"data/{prefix}_parsed_uniprot_swiss_data.csv")
            ]
            
            files_found = False
            for acc_file, sim_file, csv_file in possible_files:
                if os.path.exists(acc_file) and os.path.exists(sim_file) and os.path.exists(csv_file):
                    # Found a complete set of files
                    with st.spinner(f"Loading {dataset_type.lower()} data from local files..."):
                        df, accession_arrays, similarity_array = load_data(dataset_type)
                        if df is not None and accession_arrays is not None and similarity_array is not None:
                            st.success(f"✅ {dataset_type} data loaded from local files!")
                            files_found = True
                            break
            
            if not files_found:
                df = None
                accession_arrays = None
                similarity_array = None
    
    # Check if data loaded successfully
    if df is None or accession_arrays is None or similarity_array is None:
        st.info("💡 Please upload a zip file containing the required data files to get started.")
        return
    
    # Sidebar for controls
    st.sidebar.header("Controls")
    
    # Domain selection - use batch processing
    all_accessions = []
    for accessions in accession_arrays:
        all_accessions.extend(accessions)
    
    # Batch lookup for domains
    mask = df['accession'].isin(all_accessions)
    relevant_df = df[mask]
    
    domains = set()
    for acc in all_accessions:
        row = relevant_df[relevant_df['accession'] == acc]
        if not row.empty and pd.notna(row.iloc[0]['taxonomy_lineage']):
            taxonomy = row.iloc[0]['taxonomy_lineage']
            domain = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
            domains.add(domain)
    
    domains = sorted(list(domains))
    selected_domain = st.sidebar.selectbox(
        "Select Domain to Highlight:",
        ["All Domains"] + domains
    )
    
    # Instructions
    st.info("💡 **Tip**: Click on accession numbers in the table below to open the corresponding UniProt page in your browser!")
    
    # Main content - heatmaps on separate rows
    # Add info about synchronized sorting
    if selected_domain != "All Domains":
        st.info(f"Both heatmaps are sorted by {selected_domain} first, then by similarity within each domain.")
    
    st.subheader("Domain Heatmap")
    fig = create_heatmap(df, accession_arrays, similarity_array, selected_domain if selected_domain != "All Domains" else None)
    st.plotly_chart(fig, use_container_width=True)
    
    st.subheader("Similarity Heatmap")
    sim_fig = create_similarity_heatmap(df, accession_arrays, similarity_array, selected_domain if selected_domain != "All Domains" else None)
    st.plotly_chart(sim_fig, use_container_width=True)
    
    # Table view
    st.subheader("Protein Information Table")
    
    # Create table data with proper formatting
    table_data = []
    
    # Flatten all accessions for batch processing
    all_accessions = []
    for accessions in accession_arrays:
        all_accessions.extend(accessions)
    
    # Batch lookup using pandas isin()
    mask = df['accession'].isin(all_accessions)
    relevant_df = df[mask].set_index('accession')
    
    # If a domain is selected, sort the arrays first
    if selected_domain and selected_domain != "All Domains":
        # Sort each CDS array by domain and similarity
        sorted_arrays = []
        sorted_similarities = []
        
        for i, accessions in enumerate(accession_arrays):
            similarities = similarity_array[i].tolist() if hasattr(similarity_array[i], 'tolist') else similarity_array[i]
            
            # Get protein data for sorting using vectorized lookup
            protein_data = []
            for j, acc in enumerate(accessions):
                if acc in relevant_df.index:
                    row = relevant_df.loc[acc]
                    if pd.notna(row['taxonomy_lineage']):
                        taxonomy = row['taxonomy_lineage']
                        domain = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
                    else:
                        domain = 'Unknown'
                else:
                    domain = 'Unknown'
                protein_data.append((acc, domain, similarities[j], j))  # Include original position
            
            # Sort by domain first, then by similarity
            def sort_key(item):
                acc, domain, sim, pos = item
                if domain == selected_domain:
                    return (0, -sim)  # Selected domain first, then by similarity (descending)
                else:
                    return (1, domain, -sim)  # Other domains, sorted by domain name, then similarity
            
            sorted_data = sorted(protein_data, key=sort_key)
            sorted_accessions = [acc for acc, domain, sim, pos in sorted_data]
            sorted_sims = [sim for acc, domain, sim, pos in sorted_data]
            
            sorted_arrays.append(sorted_accessions)
            sorted_similarities.append(sorted_sims)
        
        # Create table data: iterate by position first, then by CDS (alternating CDS order)
        max_len = max(len(arr) for arr in sorted_arrays)
        for pos in range(max_len):
            for i in range(len(accession_arrays)):  # For each CDS
                sorted_accessions = sorted_arrays[i]
                if pos < len(sorted_accessions):  # Only process if this position exists in this CDS
                    acc = sorted_accessions[pos]
                    
                    # Use vectorized lookup
                    if acc in relevant_df.index:
                        row = relevant_df.loc[acc]
                        function_text = row['function']
                        if pd.notna(function_text):
                            function_display = function_text[:50] + '...' if len(function_text) > 50 else function_text
                        else:
                            function_display = 'No function data'
                        
                        # Get domain information
                        if pd.notna(row['taxonomy_lineage']):
                            taxonomy = row['taxonomy_lineage']
                            domain = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
                        else:
                            domain = 'Unknown'
                        
                        # Create clickable accession link using HTML
                        uniprot_url = f"https://www.uniprot.org/uniprot/{acc}"
                        accession_link = f'<a href="{uniprot_url}" target="_blank">{acc}</a>'
                        
                        table_data.append({
                            'CDS': f"CDS {i+1}",
                            'Accession': accession_link,
                            'Name': row['name'] if pd.notna(row['name']) else 'Not found',
                            'Full Name': row['full_name'] if pd.notna(row['full_name']) else 'Not found',
                            'Domain': domain,
                            'Organism': row['organism_common'] if pd.notna(row['organism_common']) else 'Not found',
                            'Length': row['sequence_length'] if pd.notna(row['sequence_length']) else 'Not found',
                            'Function': function_display
                        })
                    else:
                        uniprot_url = f"https://www.uniprot.org/uniprot/{acc}"
                        accession_link = f'<a href="{uniprot_url}" target="_blank">{acc}</a>'
                        table_data.append({
                            'CDS': f"CDS {i+1}",
                            'Accession': accession_link,
                            'Name': "Not found",
                            'Full Name': "Not found",
                            'Domain': "Unknown",
                            'Organism': "Not found",
                            'Length': "Not found",
                            'Function': "Not found"
                        })
    else:
        # Use original arrays without sorting
        # Get the maximum length of any CDS array
        max_len = max(len(arr) for arr in accession_arrays)
        
        # Iterate by position (index) first, then by CDS
        for pos in range(max_len):
            for i, accessions in enumerate(accession_arrays):
                if pos < len(accessions):  # Only process if this position exists in this CDS
                    acc = accessions[pos]
                    
                    # Use vectorized lookup
                    if acc in relevant_df.index:
                        row = relevant_df.loc[acc]
                        function_text = row['function']
                        if pd.notna(function_text):
                            function_display = function_text[:50] + '...' if len(function_text) > 50 else function_text
                        else:
                            function_display = 'No function data'
                        
                        # Get domain information
                        if pd.notna(row['taxonomy_lineage']):
                            taxonomy = row['taxonomy_lineage']
                            domain = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
                        else:
                            domain = 'Unknown'
                        
                        # Create clickable accession link using HTML
                        uniprot_url = f"https://www.uniprot.org/uniprot/{acc}"
                        accession_link = f'<a href="{uniprot_url}" target="_blank">{acc}</a>'
                        
                        table_data.append({
                            'CDS': f"CDS {i+1}",
                            'Accession': accession_link,
                            'Name': row['name'] if pd.notna(row['name']) else 'Not found',
                            'Full Name': row['full_name'] if pd.notna(row['full_name']) else 'Not found',
                            'Domain': domain,
                            'Organism': row['organism_common'] if pd.notna(row['organism_common']) else 'Not found',
                            'Length': row['sequence_length'] if pd.notna(row['sequence_length']) else 'Not found',
                            'Function': function_display
                        })
                    else:
                        uniprot_url = f"https://www.uniprot.org/uniprot/{acc}"
                        accession_link = f'<a href="{uniprot_url}" target="_blank">{acc}</a>'
                        table_data.append({
                            'CDS': f"CDS {i+1}",
                            'Accession': accession_link,
                            'Name': "Not found",
                            'Full Name': "Not found",
                            'Domain': "Unknown",
                            'Organism': "Not found",
                            'Length': "Not found",
                            'Function': "Not found"
                        })
    
    # Display table using HTML
    if table_data:
        # Create HTML table
        html_table = """
        <table style="width:100%; border-collapse: collapse; margin: 20px 0;">
        <thead>
        <tr style="background-color: #f0f0f0;">
        """
        
        # Add headers
        headers = ['CDS', 'Accession', 'Name', 'Full Name', 'Domain', 'Organism', 'Length', 'Function']
        for header in headers:
            html_table += f'<th style="border: 1px solid #ddd; padding: 8px; text-align: left;">{header}</th>'
        
        html_table += "</tr></thead><tbody>"
        
        # Add rows
        for row in table_data:
            html_table += '<tr>'
            for header in headers:
                cell_value = row[header]
                html_table += f'<td style="border: 1px solid #ddd; padding: 8px;">{cell_value}</td>'
            html_table += '</tr>'
        
        html_table += "</tbody></table>"
        
        # Make table scrollable
        st.markdown(f"""
        <div style="overflow-x: auto; overflow-y: auto; max-height: 500px; max-width: 100%;">
            {html_table}
        </div>
        """, unsafe_allow_html=True)
    else:
        st.write("No data to display")

if __name__ == "__main__":
    main() 