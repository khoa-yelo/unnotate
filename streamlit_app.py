"""
Streamlit version of the protein accession visualizer
"""
import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import os
import requests
import zipfile
import tempfile
from matplotlib import cm, colors as mcolors

# Helper for categorical color mapping
def get_categorical_colormap(categories):
    base_colors = list(mcolors.TABLEAU_COLORS.values()) + list(mcolors.CSS4_COLORS.values())
    color_map = {}
    for i, cat in enumerate(categories):
        color_map[cat] = base_colors[i % len(base_colors)]
    return color_map

def download_from_github_release(repo_owner, repo_name, release_tag="latest"):
    """Download data files from GitHub release if they don't exist locally"""
    
    # Check if files already exist
    required_files = [
        "full_accession_arrays.npy",
        "full_similarity_array.npy", 
        "parsed_uniprot_swiss_data.csv"
    ]
    
    existing_files = [f for f in required_files if os.path.exists(f)]
    if len(existing_files) == len(required_files):
        st.success("✅ All data files found locally")
        return True
    
    st.info(f"📥 Downloading missing data files from GitHub release...")
    
    if release_tag == "latest":
        api_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/releases/latest"
    else:
        api_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/releases/tags/{release_tag}"
    
    try:
        response = requests.get(api_url)
        response.raise_for_status()
        release_data = response.json()
        
        # Find the zip file asset
        zip_asset = None
        for asset in release_data.get('assets', []):
            if asset['name'].endswith('.zip'):
                zip_asset = asset
                break
        
        if not zip_asset:
            st.error("No zip file found in the release")
            return False
        
        # Download the zip file
        zip_url = zip_asset['browser_download_url']
        
        with st.spinner("Downloading data files..."):
            zip_response = requests.get(zip_url, stream=True)
            zip_response.raise_for_status()
            
            # Save to temporary file
            with tempfile.NamedTemporaryFile(delete=False, suffix='.zip') as tmp_file:
                for chunk in zip_response.iter_content(chunk_size=8192):
                    tmp_file.write(chunk)
                tmp_path = tmp_file.name
            
            # Extract to current directory
            with zipfile.ZipFile(tmp_path, 'r') as zip_ref:
                zip_ref.extractall('.')
            
            # Clean up temporary file
            os.unlink(tmp_path)
        
        st.success("✅ Data files downloaded and extracted successfully!")
        return True
        
    except Exception as e:
        st.error(f"Failed to download from GitHub: {e}")
        return False

@st.cache_data
def load_data():
    """Load and cache the protein data"""
    # Try multiple possible paths for the CSV file
    csv_paths = [
        "parsed_uniprot_swiss_data.csv",
        "data/unnotate/parsed_uniprot_swiss_data.csv",
        "../data/unnotate/parsed_uniprot_swiss_data.csv"
    ]
    
    df = None
    for path in csv_paths:
        if os.path.exists(path):
            df = pd.read_csv(path)
            break
    
    if df is None:
        st.error("❌ Could not find parsed_uniprot_swiss_data.csv")
        st.info("Please download the data files using the data downloader or upload them manually.")
        return None, None, None
    
    # Try to load numpy arrays
    accession_arrays = None
    similarity_array = None
    
    # Try multiple possible paths for numpy files
    npy_paths = [
        ("full_accession_arrays.npy", "full_similarity_array.npy"),
        ("accession_arrays.npy", "similarity_array.npy"),
        ("viral_accession_arrays.npy", "viral_similarity_array.npy"),
        ("bacterial_accession_arrays.npy", "bacterial_similarity_array.npy")
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
                
                st.success(f"✅ Loaded data from {acc_path} and {sim_path}")
                break
                
            except Exception as e:
                st.warning(f"Failed to load {acc_path}: {e}")
                continue
    
    if accession_arrays is None:
        st.warning("⚠️ Could not load numpy arrays, using sample data")
        # Fallback to sample data
        accession_arrays = [
            ['P03756', 'Q21LK2', 'Q89EM9', 'Q03544', 'Q03546', 'P76515', 'A4IM80', 'B6JPK6', 'B5Z9N6', 'Q03547'],
            ['P43661', 'Q887Q8', 'Q8X5K6', 'P59791', 'P33128', 'Q88ND4', 'P70799', 'Q06062', 'P39632', 'P33409']
        ]
        similarity_array = np.array([
            [1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.3, 0.5, 0.7, 0.9],
            [0.8, 1.0, 0.7, 0.5, 0.3, 0.2, 0.4, 0.6, 0.8, 0.7]
        ])
    
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
    
    # Check if data files exist, if not try to download from GitHub
    required_files = [
        "full_accession_arrays.npy",
        "full_similarity_array.npy", 
        "parsed_uniprot_swiss_data.csv"
    ]
    
    missing_files = [f for f in required_files if not os.path.exists(f)]
    
    if missing_files:
        st.warning(f"⚠️ Missing data files: {', '.join(missing_files)}")
        
        # Auto-download for cloud deployment
        if st.secrets.get("GITHUB_REPO_OWNER") and st.secrets.get("GITHUB_REPO_NAME"):
            # Use secrets for cloud deployment
            repo_owner = st.secrets["GITHUB_REPO_OWNER"]
            repo_name = st.secrets["GITHUB_REPO_NAME"]
            release_tag = st.secrets.get("GITHUB_RELEASE_TAG", "latest")
            
            st.info("🔄 Auto-downloading data from GitHub release...")
            success = download_from_github_release(repo_owner, repo_name, release_tag)
            if success:
                st.rerun()
        else:
            # Manual download interface for local development
            with st.expander("📥 Download Data from GitHub Release", expanded=True):
                st.write("The app needs data files to run. You can download them from the GitHub release:")
                
                col1, col2 = st.columns(2)
                with col1:
                    repo_owner = st.text_input("Repository Owner", value="khoa-yelo")
                    repo_name = st.text_input("Repository Name", value="unnotate")
                
                with col2:
                    release_tag = st.text_input("Release Tag", value="latest", 
                                               help="Use 'latest' for the most recent release, or specify a tag like 'v1.0.1'")
                    
                    if st.button("Download Data"):
                        if repo_owner and repo_name:
                            success = download_from_github_release(repo_owner, repo_name, release_tag)
                            if success:
                                st.rerun()  # Refresh the page
                        else:
                            st.error("Please enter repository owner and name")
                
                st.info("💡 **Alternative**: You can also manually download the `protein_data_v1.0.zip` file from the GitHub releases page and extract it here.")
    
    # Load data
    with st.spinner("Loading data..."):
        df, accession_arrays, similarity_array = load_data()
    
    # Check if data loaded successfully
    if df is None or accession_arrays is None or similarity_array is None:
        st.error("❌ Failed to load required data.")
        st.info("💡 Please download the data files using the downloader above or upload them manually.")
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