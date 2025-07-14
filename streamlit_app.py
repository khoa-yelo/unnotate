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

def precompute_domain_data(df, accession_arrays):
    """Precompute domain information and hover text for all accessions"""
    # Flatten all accessions
    all_accessions = []
    for accessions in accession_arrays:
        all_accessions.extend(accessions)
    
    # Batch lookup using pandas isin() - do this once
    mask = df['accession'].isin(all_accessions)
    relevant_df = df[mask].set_index('accession')
    
    # Precompute domain and hover data for each accession
    domain_cache = {}
    hover_cache = {}
    
    for acc in all_accessions:
        if acc in relevant_df.index:
            row = relevant_df.loc[acc]
            if pd.notna(row['taxonomy_lineage']):
                taxonomy = row['taxonomy_lineage']
                domain = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
            else:
                domain = 'Unknown'
            
            # Precompute hover text
            function_text = row['function']
            if pd.notna(function_text):
                function_display = function_text[:50] + ('...' if len(function_text) > 50 else '')
            else:
                function_display = 'No function data'
            
            hover_text = f"""
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
        
        domain_cache[acc] = domain
        hover_cache[acc] = hover_text
    
    return domain_cache, hover_cache, relevant_df

def precompute_sorted_arrays(accession_arrays, similarity_array, domain_cache, selected_domain):
    """Precompute sorted arrays for a given domain selection"""
    if not selected_domain or selected_domain == "All Domains":
        return accession_arrays, similarity_array
    
    sorted_arrays = []
    sorted_similarities = []
    
    for i, accessions in enumerate(accession_arrays):
        similarities = similarity_array[i].tolist() if hasattr(similarity_array[i], 'tolist') else similarity_array[i]
        
        # Get protein data for sorting using cached domain data
        protein_data = []
        for j, acc in enumerate(accessions):
            domain = domain_cache.get(acc, 'Unknown')
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
    
    return sorted_arrays, np.array(sorted_similarities)

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
                sequence_similarity_file = None
                
                for file_path in file_list:
                    if file_path.endswith('.csv'):
                        csv_file = os.path.join(temp_dir, file_path)
                    elif 'accession' in file_path.lower() and file_path.endswith('.npy'):
                        accession_file = os.path.join(temp_dir, file_path)
                    elif 'similarity' in file_path.lower() and 'sequence' not in file_path.lower() and file_path.endswith('.npy'):
                        similarity_file = os.path.join(temp_dir, file_path)
                    elif 'sequence_similarity' in file_path.lower() and file_path.endswith('.npy'):
                        sequence_similarity_file = os.path.join(temp_dir, file_path)
                
                if not csv_file:
                    st.error("❌ No CSV file found in the zip")
                    return None, None, None, None
                
                if not accession_file:
                    st.error("❌ No accession arrays file found in the zip")
                    return None, None, None, None
                
                if not similarity_file:
                    st.error("❌ No similarity array file found in the zip")
                    return None, None, None, None
                
                # Load the data
                df = pd.read_csv(csv_file)
                st.success(f"✅ Loaded CSV with {len(df)} rows")
                
                accession_arrays_loaded = np.load(accession_file)
                similarity_array = np.load(similarity_file)
                
                # Load sequence similarity array if available
                sequence_similarity_array = None
                if sequence_similarity_file:
                    sequence_similarity_array = np.load(sequence_similarity_file)
                    st.success(f"✅ Loaded sequence similarity array with shape {sequence_similarity_array.shape}")
                
                # Convert to list format
                accession_arrays = []
                for arr in accession_arrays_loaded:
                    if hasattr(arr, 'tolist'):
                        arr = arr.tolist()
                    arr = [acc.decode('utf-8') if isinstance(acc, bytes) else acc for acc in arr]
                    accession_arrays.append(arr)
                
                st.success(f"✅ Loaded {len(accession_arrays)} accession arrays and similarity array")
                
                return df, accession_arrays, similarity_array, sequence_similarity_array
                
    except Exception as e:
        st.error(f"❌ Error processing zip file: {e}")
        return None, None, None, None

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
        return None, None, None, None
    
    # Try to load numpy arrays
    accession_arrays = None
    similarity_array = None
    sequence_similarity_array = None
    
    # Try multiple possible paths for numpy files
    npy_paths = [
        (f"{prefix}_accession_arrays.npy", f"{prefix}_similarity_array.npy", f"{prefix}_sequence_similarity_array.npy"),
        (f"data/{prefix}_accession_arrays.npy", f"data/{prefix}_similarity_array.npy", f"data/{prefix}_sequence_similarity_array.npy"),
        (str(data_path / f"{prefix}_accession_arrays.npy"), str(data_path / f"{prefix}_similarity_array.npy"), str(data_path / f"{prefix}_sequence_similarity_array.npy")),
        (f"../data/{prefix}_accession_arrays.npy", f"../data/{prefix}_similarity_array.npy", f"../data/{prefix}_sequence_similarity_array.npy"),
    ]
    
    for acc_path, sim_path, seq_sim_path in npy_paths:
        if os.path.exists(acc_path) and os.path.exists(sim_path):
            try:
                accession_arrays_loaded = np.load(acc_path)
                similarity_array = np.load(sim_path)
                
                # Load sequence similarity array if available
                if os.path.exists(seq_sim_path):
                    sequence_similarity_array = np.load(seq_sim_path)
                
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
        return None, None, None, None
    
    return df, accession_arrays, similarity_array, sequence_similarity_array

def create_heatmap(df, accession_arrays, similarity_array, selected_domain=None, domain_cache=None, hover_cache=None):
    """Create heatmap visualization with optional domain-based sorting"""
    # Use precomputed data if available, otherwise compute it
    if domain_cache is None or hover_cache is None:
        domain_cache, hover_cache, relevant_df = precompute_domain_data(df, accession_arrays)
    
    # Get sorted arrays using precomputed function
    arrays_to_visualize, sorted_similarity_array = precompute_sorted_arrays(
        accession_arrays, similarity_array, domain_cache, selected_domain
    )
    
    max_len = max(len(arr) for arr in arrays_to_visualize)
    
    # Create domain data efficiently using precomputed caches
    domain_data = []
    hover_texts = []
    customdata = []
    
    for cds_idx, accessions in enumerate(arrays_to_visualize):
        row_domains = []
        row_hover = []
        row_urls = []
        
        for pos_idx, acc in enumerate(accessions):
            # Use precomputed domain and hover data
            domain = domain_cache.get(acc, 'Unknown')
            base_hover = hover_cache.get(acc, f"<b>Accession:</b> {acc}<br><b>Status:</b> Not found in database<br><b>Click to open UniProt page</b>")
            
            # Add CDS and position info to hover text
            hover_text = f"""
            <b>CDS Index:</b> {cds_idx+1}<br>
            <b>Position Index:</b> {pos_idx+1}<br>
            {base_hover}
            """
            
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
        yaxis_title="Protein Index",
        height=400,
        width=1200,
        legend_title="Domain"
    )
    return fig

def create_similarity_heatmap(df, accession_arrays, similarity_array, selected_domain=None, domain_cache=None, hover_cache=None):
    """Create similarity heatmap with optional domain-based sorting"""
    # Use precomputed data if available, otherwise compute it
    if domain_cache is None or hover_cache is None:
        domain_cache, hover_cache, _ = precompute_domain_data(df, accession_arrays)
    
    # Get sorted arrays using precomputed function
    arrays_to_visualize, similarity_data = precompute_sorted_arrays(
        accession_arrays, similarity_array, domain_cache, selected_domain
    )
    
    # Create hover text and customdata with UniProt URLs
    max_len = max(len(arr) for arr in arrays_to_visualize)
    hover_texts = []
    customdata = []
    
    for cds_idx, accessions in enumerate(arrays_to_visualize):
        row_hover = []
        row_urls = []
        
        for pos_idx, acc in enumerate(accessions):
            # Get similarity value
            similarity_value = similarity_data[cds_idx][pos_idx]
            
            # Use precomputed hover data
            base_hover = hover_cache.get(acc, f"<b>Accession:</b> {acc}<br><b>Status:</b> Not found in database<br><b>Click to open UniProt page</b>")
            
            # Add CDS, position, and similarity info to hover text
            hover_text = f"""
            <b>CDS Index:</b> {cds_idx+1}<br>
            <b>Position Index:</b> {pos_idx+1}<br>
            <b>Similarity:</b> {similarity_value:.3f}<br>
            {base_hover}
            """
            
            uniprot_url = f"https://www.uniprot.org/uniprot/{acc}"
            
            row_hover.append(hover_text)
            row_urls.append(uniprot_url)
        
        # Pad with empty strings
        row_hover.extend([''] * (max_len - len(accessions)))
        row_urls.extend([''] * (max_len - len(accessions)))
        hover_texts.append(row_hover)
        customdata.append(row_urls)
    
    # Transpose and flip
    hover_texts_transposed = list(zip(*hover_texts))[::-1]
    customdata_transposed = list(zip(*customdata))[::-1]
    
    z = np.array(similarity_data).T[::-1]
    
    # Normalize values if they're in 0-100 range
    if np.max(z) > 1.1:
        z = z / 100.0
    
    fig = go.Figure(data=go.Heatmap(
        z=z,
        colorscale='Blues',
        colorbar=dict(title="Similarity"),
        hovertext=hover_texts_transposed,
        hoverinfo='text',
        customdata=customdata_transposed,
        zmin=0.8,
        zmax=1
    ))
    
    title = "Cosine Similarity (ESM-C) Heatmap"
    if selected_domain and selected_domain != "All Domains":
        title += f" - Sorted by {selected_domain} first, then by similarity"
    
    fig.update_layout(
        title=title,
        xaxis_title="CDS Index",
        yaxis_title="Protein Index",    
        height=400,
        width=1200
    )
    return fig

def create_similarity_barplot(df, accession_arrays, similarity_array, selected_domain=None, domain_cache=None):
    """Create horizontal bar plot showing row sums for similarity data"""
    # Use precomputed data if available, otherwise compute it
    if domain_cache is None:
        domain_cache, _, _ = precompute_domain_data(df, accession_arrays)
    
    # Get sorted arrays using precomputed function
    arrays_to_visualize, similarity_data = precompute_sorted_arrays(
        accession_arrays, similarity_array, domain_cache, selected_domain
    )
    
    # Calculate row sums
    z = np.array(similarity_data).T[::-1]
    
    # Normalize values if they're in 0-100 range
    if np.max(z) > 1.1:
        z = z / 100.0
    
    row_sums = np.sum(z, axis=1)
    
    # Create bar plot
    fig = go.Figure(data=go.Bar(
        x=row_sums,
        y=[f"Row {i+1}" for i in range(len(row_sums))],
        orientation='h',
        marker_color='lightblue',
        marker_line_color='navy',
        marker_line_width=1,
        text=[f"{val:.2f}" for val in row_sums],
        textposition='auto',
        hovertemplate='<b>Row %{y}</b><br><b>Total Similarity:</b> %{x:.3f}<br><b>Average:</b> %{customdata:.3f}<extra></extra>',
        customdata=[val/len(z[0]) for val in row_sums]
    ))
    
    title = "Row Sums - Cosine Similarity"
    if selected_domain and selected_domain != "All Domains":
        title += f" (Sorted by {selected_domain})"
    
    fig.update_layout(
        title=title,
        xaxis_title="Total Similarity",
        yaxis_title="Protein Row",
        height=400,
        width=300,
        showlegend=False
    )
    
    return fig

def create_sequence_similarity_heatmap(df, accession_arrays, sequence_similarity_array, selected_domain=None, domain_cache=None, hover_cache=None):
    """Create Sequence  Identity Heatmap with optional domain-based sorting"""
    # Use precomputed data if available, otherwise compute it
    if domain_cache is None or hover_cache is None:
        domain_cache, hover_cache, _ = precompute_domain_data(df, accession_arrays)
    
    # Get sorted arrays using precomputed function
    arrays_to_visualize, sequence_similarity_data = precompute_sorted_arrays(
        accession_arrays, sequence_similarity_array, domain_cache, selected_domain
    )
    
    # Create hover text and customdata with UniProt URLs
    max_len = max(len(arr) for arr in arrays_to_visualize)
    hover_texts = []
    customdata = []
    
    for cds_idx, accessions in enumerate(arrays_to_visualize):
        row_hover = []
        row_urls = []
        
        for pos_idx, acc in enumerate(accessions):
            # Get sequence similarity value
            seq_similarity_value = sequence_similarity_data[cds_idx][pos_idx]
            
            # Use precomputed hover data
            base_hover = hover_cache.get(acc, f"<b>Accession:</b> {acc}<br><b>Status:</b> Not found in database<br><b>Click to open UniProt page</b>")
            
            # Add CDS, position, and sequence similarity info to hover text
            hover_text = f"""
            <b>CDS Index:</b> {cds_idx+1}<br>
            <b>Position Index:</b> {pos_idx+1}<br>
            <b>Sequence Identity:</b> {seq_similarity_value:.3f}<br>
            {base_hover}
            """
            
            uniprot_url = f"https://www.uniprot.org/uniprot/{acc}"
            
            row_hover.append(hover_text)
            row_urls.append(uniprot_url)
        
        # Pad with empty strings
        row_hover.extend([''] * (max_len - len(accessions)))
        row_urls.extend([''] * (max_len - len(accessions)))
        hover_texts.append(row_hover)
        customdata.append(row_urls)
    
    # Transpose and flip
    hover_texts_transposed = list(zip(*hover_texts))[::-1]
    customdata_transposed = list(zip(*customdata))[::-1]
    
    z = np.array(sequence_similarity_data).T[::-1]
    
    # Normalize values if they're in 0-100 range
    if np.max(z) > 1.1:
        z = z / 100.0
    
    fig = go.Figure(data=go.Heatmap(
        z=z,
        colorscale='Reds',
        colorbar=dict(title="Sequence Similarity"),
        hovertext=hover_texts_transposed,
        hoverinfo='text',
        customdata=customdata_transposed,
        zmin=0,
        zmax=1
    ))
    
    title = "Sequence  Identity Heatmap"
    if selected_domain and selected_domain != "All Domains":
        title += f" - Sorted by {selected_domain} first, then by similarity"
    
    fig.update_layout(
        title=title,
        xaxis_title="CDS Index",
        yaxis_title="Protein Index",    
        height=400,
        width=1200
    )
    return fig

def create_sequence_barplot(df, accession_arrays, sequence_similarity_array, selected_domain=None, domain_cache=None):
    """Create horizontal bar plot showing row sums for sequence similarity data"""
    # Use precomputed data if available, otherwise compute it
    if domain_cache is None:
        domain_cache, _, _ = precompute_domain_data(df, accession_arrays)
    
    # Get sorted arrays using precomputed function
    arrays_to_visualize, sequence_similarity_data = precompute_sorted_arrays(
        accession_arrays, sequence_similarity_array, domain_cache, selected_domain
    )
    
    # Calculate row sums
    z = np.array(sequence_similarity_data).T[::-1]
    
    # Normalize values if they're in 0-100 range
    if np.max(z) > 1.1:
        z = z / 100.0
    
    row_sums = np.sum(z, axis=1)
    
    # Create bar plot
    fig = go.Figure(data=go.Bar(
        x=row_sums,
        y=[f"Row {i+1}" for i in range(len(row_sums))],
        orientation='h',
        marker_color='lightcoral',
        marker_line_color='darkred',
        marker_line_width=1,
        text=[f"{val:.2f}" for val in row_sums],
        textposition='auto',
        hovertemplate='<b>Row %{y}</b><br><b>Total Sequence Identity:</b> %{x:.3f}<br><b>Average:</b> %{customdata:.3f}<extra></extra>',
        customdata=[val/len(z[0]) for val in row_sums]
    ))
    
    title = "Row Sums - Sequence Identity"
    if selected_domain and selected_domain != "All Domains":
        title += f" (Sorted by {selected_domain})"
    
    fig.update_layout(
        title=title,
        xaxis_title="Total Sequence Identity",
        yaxis_title="Protein Row",
        height=400,
        width=300,
        showlegend=False
    )
    
    return fig

def create_similarity_scatter_plot(df, accession_arrays, similarity_array, sequence_similarity_array, selected_domain=None, domain_cache=None, hover_cache=None):
    """Create interactive scatter plot comparing cosine similarity vs sequence identity"""
    # Use precomputed data if available, otherwise compute it
    if domain_cache is None or hover_cache is None:
        domain_cache, hover_cache, _ = precompute_domain_data(df, accession_arrays)
    
    # Prepare data for scatter plot
    cosine_similarities = []
    sequence_identities = []
    accessions = []
    hover_texts = []
    colors = []
    
    # Get domain colors for consistent coloring
    unique_domains = set()
    for accessions_list in accession_arrays:
        for acc in accessions_list:
            domain = domain_cache.get(acc, 'Unknown')
            unique_domains.add(domain)
    
    domain_colors = get_categorical_colormap(unique_domains)
    
    for cds_idx, accessions_list in enumerate(accession_arrays):
        for pos_idx, acc in enumerate(accessions_list):
            # Get similarity values
            cosine_sim = similarity_array[cds_idx][pos_idx]
            seq_id = sequence_similarity_array[cds_idx][pos_idx]
            
            # Normalize values if they're in 0-100 range
            if cosine_sim > 1.1:
                cosine_sim = cosine_sim / 100.0
            if seq_id > 1.1:
                seq_id = seq_id / 100.0
            
            # Get domain for coloring
            domain = domain_cache.get(acc, 'Unknown')
            color = domain_colors.get(domain, '#808080')  # Default gray
            
            # Create hover text
            base_hover = hover_cache.get(acc, f"<b>Accession:</b> {acc}<br><b>Status:</b> Not found in database<br><b>Click to open UniProt page</b>")
            hover_text = f"""
            <b>CDS Index:</b> {cds_idx+1}<br>
            <b>Position Index:</b> {pos_idx+1}<br>
            <b>Cosine Similarity:</b> {cosine_sim:.3f}<br>
            <b>Sequence Identity:</b> {seq_id:.3f}<br>
            <b>Domain:</b> {domain}<br>
            {base_hover}
            """
            
            cosine_similarities.append(cosine_sim)
            sequence_identities.append(seq_id)
            accessions.append(acc)
            hover_texts.append(hover_text)
            colors.append(color)
    
    # Create scatter plot
    fig = go.Figure()
    
    # Add scatter trace
    fig.add_trace(go.Scatter(
        x=cosine_similarities,
        y=sequence_identities,
        mode='markers',
        marker=dict(
            size=8,
            color=colors,
            opacity=0.7,
            line=dict(width=1, color='white')
        ),
        text=hover_texts,
        hoverinfo='text',
        customdata=[f"https://www.uniprot.org/uniprot/{acc}" for acc in accessions],
        hovertemplate='%{text}<extra></extra>'
    ))
    
    # Calculate correlation coefficient
    correlation = np.corrcoef(cosine_similarities, sequence_identities)[0, 1]
    
    # Update layout
    fig.update_layout(
        title=f"Cosine Similarity vs Sequence Identity (Correlation: {correlation:.3f})",
        xaxis_title="Cosine Similarity (ESM-C)",
        yaxis_title="Sequence Identity",
        xaxis=dict(range=[0.8, 1.0]),  # Match cosine similarity heatmap zmin/zmax
        yaxis=dict(range=[0, 1]),       # Match sequence identity heatmap zmin/zmax
        height=500,
        width=800,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        )
    )
    
    # Add domain legend
    for domain in sorted(unique_domains):
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(size=10, color=domain_colors.get(domain, '#808080')),
            name=domain,
            showlegend=True
        ))
    
    return fig

def main():
    st.set_page_config(page_title="Protein Accession Visualizer", layout="wide")
    
    st.title("Unnotate: Annotation of Proteins with Unknown function + Uncertainty quantification + Uniprot Mapping")
    
    # Dataset selection dropdown
    st.sidebar.header("Dataset Selection")
    dataset_type = st.sidebar.selectbox(
        "Choose Dataset:",
        ["Virus", "Bacteria"],
        help="Select which default  dataset to visualize"
    )
    
    # Always show upload interface
    with st.expander("📥 Upload Data from ZIP File", expanded=True):
        st.write("Upload a zip file containing the required data files:")
        st.write("- A CSV file (e.g., `PREFIX_parsed_uniprot_swiss_data.csv`)")
        st.write("- A numpy file with accession arrays (e.g., `PREFIX_accession_arrays.npy`)")
        st.write("- A numpy file with similarity data (e.g., `PREFIX_similarity_array.npy`)")
        st.write("- A numpy file with sequence similarity data (e.g., `PREFIX_sequence_similarity_array.npy`) - optional")
        
        uploaded_file = st.file_uploader("Choose a ZIP file", type=["zip"])
        
        if uploaded_file is not None:
            with st.spinner("Processing uploaded ZIP file..."):
                df, accession_arrays, similarity_array, sequence_similarity_array = process_uploaded_zip(uploaded_file)
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
                        df, accession_arrays, similarity_array, sequence_similarity_array = load_data(dataset_type)
                        if df is not None and accession_arrays is not None and similarity_array is not None:
                            st.success(f"✅ {dataset_type} data loaded from local files!")
                            files_found = True
                            break
            
            if not files_found:
                df = None
                accession_arrays = None
                similarity_array = None
                sequence_similarity_array = None
    
    # Check if data loaded successfully
    if df is None or accession_arrays is None or similarity_array is None:
        st.info("💡 Please upload a zip file containing the required data files to get started.")
        return
    
    # Initialize session state for caching
    if 'domain_cache' not in st.session_state:
        st.session_state.domain_cache = None
        st.session_state.hover_cache = None
        st.session_state.relevant_df = None
        st.session_state.domains = None
    
    # Precompute domain data if not cached or if data changed
    data_hash = hash(str(accession_arrays) + str(df.shape))
    if (st.session_state.domain_cache is None or 
        'data_hash' not in st.session_state or 
        st.session_state.data_hash != data_hash):
        
        with st.spinner("Precomputing domain data..."):
            st.session_state.domain_cache, st.session_state.hover_cache, st.session_state.relevant_df = precompute_domain_data(df, accession_arrays)
            st.session_state.data_hash = data_hash
            
            # Precompute domains list
            all_accessions = []
            for accessions in accession_arrays:
                all_accessions.extend(accessions)
            
            mask = st.session_state.relevant_df.index.isin(all_accessions)
            relevant_subset = st.session_state.relevant_df[mask]
            
            domains = set()
            for acc in all_accessions:
                if acc in relevant_subset.index:
                    row = relevant_subset.loc[acc]
                    if pd.notna(row['taxonomy_lineage']):
                        taxonomy = row['taxonomy_lineage']
                        domain = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
                        domains.add(domain)
            
            st.session_state.domains = sorted(list(domains))
    
    # Sidebar for controls
    st.sidebar.header("Controls")
    
    # Domain selection using cached domains
    domains = st.session_state.domains or ["All Domains"]
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
    fig = create_heatmap(
        df, accession_arrays, similarity_array, 
        selected_domain if selected_domain != "All Domains" else None,
        st.session_state.domain_cache, st.session_state.hover_cache
    )
    st.plotly_chart(fig, use_container_width=True)
    
    # Cosine Similarity Heatmap and Bar Plot
    st.subheader("Cosine Similarity Heatmap (ESM-C)")
    col1, col2 = st.columns([4, 1])
    
    with col1:
        sim_fig = create_similarity_heatmap(
            df, accession_arrays, similarity_array, 
            selected_domain if selected_domain != "All Domains" else None,
            st.session_state.domain_cache,
            st.session_state.hover_cache
        )
        st.plotly_chart(sim_fig, use_container_width=True)
    
    with col2:
        sim_bar_fig = create_similarity_barplot(
            df, accession_arrays, similarity_array,
            selected_domain if selected_domain != "All Domains" else None,
            st.session_state.domain_cache
        )
        st.plotly_chart(sim_bar_fig, use_container_width=True)
    
    # Display Sequence Identity Heatmap and Bar Plot if available
    if sequence_similarity_array is not None:
        st.subheader("Sequence Identity Heatmap")
        col3, col4 = st.columns([4, 1])
        
        with col3:
            seq_sim_fig = create_sequence_similarity_heatmap(
                df, accession_arrays, sequence_similarity_array, 
                selected_domain if selected_domain != "All Domains" else None,
                st.session_state.domain_cache,
                st.session_state.hover_cache
            )
            st.plotly_chart(seq_sim_fig, use_container_width=True)
        
        with col4:
            seq_bar_fig = create_sequence_barplot(
                df, accession_arrays, sequence_similarity_array,
                selected_domain if selected_domain != "All Domains" else None,
                st.session_state.domain_cache
            )
            st.plotly_chart(seq_bar_fig, use_container_width=True)
        
        # Add scatter plot comparing the two similarity measures
        st.subheader("Similarity Comparison Scatter Plot")
        scatter_fig = create_similarity_scatter_plot(
            df, accession_arrays, similarity_array, sequence_similarity_array,
            selected_domain if selected_domain != "All Domains" else None,
            st.session_state.domain_cache,
            st.session_state.hover_cache
        )
        st.plotly_chart(scatter_fig, use_container_width=True)
    else:
        st.info("💡 No sequence similarity data available. Upload a zip file containing 'sequence_similarity_array.npy' to see this visualization.")
    
    # Table view
    st.subheader("Protein Information Table")
    
    # Create table data with proper formatting using cached data
    table_data = []
    
    # Get sorted arrays using precomputed function
    arrays_to_visualize, sorted_similarity_array = precompute_sorted_arrays(
        accession_arrays, similarity_array, st.session_state.domain_cache, selected_domain
    )
    
    # Use cached relevant_df
    relevant_df = st.session_state.relevant_df
    
    # Create table data: iterate by position first, then by CDS (alternating CDS order)
    max_len = max(len(arr) for arr in arrays_to_visualize)
    for pos in range(max_len):
        for i in range(len(arrays_to_visualize)):  # For each CDS
            accessions = arrays_to_visualize[i]
            if pos < len(accessions):  # Only process if this position exists in this CDS
                acc = accessions[pos]
                
                # Use cached data
                if relevant_df is not None and acc in relevant_df.index:
                    row = relevant_df.loc[acc]
                    function_text = row['function']
                    if pd.notna(function_text):
                        function_display = function_text[:50] + '...' if len(function_text) > 50 else function_text
                    else:
                        function_display = 'No function data'
                    
                    # Get domain information from cache
                    domain = st.session_state.domain_cache.get(acc, 'Unknown')
                    
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