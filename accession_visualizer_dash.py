import dash
from dash import dcc, html, dash_table, Input, Output
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import os
import h5py
from protein_embedder import embed_proteins
from faiss_knn import FaissKNN

# Helper for categorical color mapping
from matplotlib import cm, colors as mcolors

def get_categorical_colormap(categories):
    base_colors = list(mcolors.TABLEAU_COLORS.values()) + list(mcolors.CSS4_COLORS.values())
    color_map = {}
    for i, cat in enumerate(categories):
        color_map[cat] = base_colors[i % len(base_colors)]
    return color_map

def get_available_kingdoms(df, accession_arrays):
    kingdoms = set()
    for accessions in accession_arrays:
        accessions = [acc.decode('utf-8') if isinstance(acc, bytes) else acc for acc in accessions]
        for acc in accessions:
            row = df[df['accession'] == acc]
            if not row.empty and pd.notna(row.iloc[0]['taxonomy_lineage']):
                taxonomy = row.iloc[0]['taxonomy_lineage']
                kingdom = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
                kingdoms.add(kingdom)
    return sorted(list(kingdoms))

def create_accession_heatmap(df, accession_arrays):
    max_len = max(len(arr) for arr in accession_arrays)
    accession_labels = []
    kingdom_data = []
    hover_texts = []
    for i, accessions in enumerate(accession_arrays):
        accessions = [acc.decode('utf-8') if isinstance(acc, bytes) else acc for acc in accessions]
        padded_accessions = accessions + [''] * (max_len - len(accessions))
        accession_labels.append(padded_accessions)
        row_kingdoms = []
        row_hover = []
        for acc in accessions:
            row = df[df['accession'] == acc]
            if not row.empty and pd.notna(row.iloc[0]['taxonomy_lineage']):
                taxonomy = row.iloc[0]['taxonomy_lineage']
                kingdom = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
            else:
                kingdom = 'Unknown'
            row_kingdoms.append(kingdom)
            if not row.empty:
                function_text = row.iloc[0]['function']
                if pd.notna(function_text):
                    function_display = function_text[:50] + ('...' if len(function_text) > 50 else '')
                else:
                    function_display = 'No function data'
                taxonomy = row.iloc[0]['taxonomy_lineage'] if pd.notna(row.iloc[0]['taxonomy_lineage']) else 'Unknown'
                kingdom = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
                full_name = row.iloc[0]['full_name'] if pd.notna(row.iloc[0]['full_name']) else 'Not found'
                hover_text = f"""
                <b>Accession:</b> {acc}<br>
                <b>Full Name:</b> {full_name}<br>
                <b>Name:</b> {row.iloc[0]['name']}<br>
                <b>Organism:</b> {row.iloc[0]['organism_common']}<br>
                <b>Kingdom:</b> {kingdom}<br>
                <b>Function:</b> {function_display}
                """
            else:
                hover_text = f"<b>Accession:</b> {acc}<br><b>Status:</b> Not found in database"
            row_hover.append(hover_text)
        row_kingdoms.extend([''] * (max_len - len(accessions)))
        row_hover.extend([''] * (max_len - len(accessions)))
        kingdom_data.append(row_kingdoms)
        hover_texts.append(row_hover)
    accession_labels_transposed = list(zip(*accession_labels))[::-1]
    hover_texts_transposed = list(zip(*hover_texts))[::-1]
    kingdom_data_transposed = list(zip(*kingdom_data))[::-1]
    unique_kingdoms = sorted(set(k for row in kingdom_data_transposed for k in row if k and k != ''))
    kingdom_to_num = {kingdom: i for i, kingdom in enumerate(unique_kingdoms)}
    kingdom_numerical = []
    for row in kingdom_data_transposed:
        num_row = []
        for kingdom in row:
            if kingdom and kingdom != '':
                num_row.append(kingdom_to_num[kingdom])
            else:
                num_row.append(-1)
        kingdom_numerical.append(num_row)
    kingdom_color_map = get_categorical_colormap(unique_kingdoms)
    colorscale = []
    for i, k in enumerate(unique_kingdoms):
        color = kingdom_color_map[k]
        if len(unique_kingdoms) == 1:
            colorscale = [[0, color], [1, color]]
        else:
            val = i / (len(unique_kingdoms) - 1)
            colorscale.append([val, color])
    if len(unique_kingdoms) > 1:
        last_color = kingdom_color_map[unique_kingdoms[-1]]
        colorscale.append([1.0, last_color])
    fig = go.Figure(data=go.Heatmap(
        z=kingdom_numerical,
        hovertext=hover_texts_transposed,
        hoverinfo='text',
        colorscale=colorscale,
        showscale=False,
        zmin=0,
        zmax=len(unique_kingdoms)-1
    ))
    for i, k in enumerate(unique_kingdoms):
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(size=15, color=kingdom_color_map[k]),
            legendgroup=k,
            showlegend=True,
            name=k
        ))
    fig.update_layout(
        title="Kingdom Heatmap",
        xaxis_title="Array Index",
        yaxis_title="Position in Array",
        height=400,
        width=1200,
        legend_title="Kingdom"
    )
    return fig

def create_similarity_heatmap(similarity_array):
    z = np.array(similarity_array).T[::-1]
    fig = go.Figure(data=go.Heatmap(
        z=z,
        colorscale='Blues',
        colorbar=dict(title="Similarity")
    ))
    fig.update_layout(
        title="Similarity Heatmap",
        xaxis_title="Array Index",
        yaxis_title="Position in Array",
        height=400,
        width=1200
    )
    return fig

def create_accession_table(df, accession_arrays):
    table_data = []
    for i, accessions in enumerate(accession_arrays):
        accessions = [acc.decode('utf-8') if isinstance(acc, bytes) else acc for acc in accessions]
        for j, acc in enumerate(accessions):
            row = df[df['accession'] == acc]
            if not row.empty:
                function_text = row.iloc[0]['function']
                if pd.notna(function_text):
                    function_display = function_text[:50] + '...' if len(function_text) > 50 else function_text
                else:
                    function_display = 'No function data'
                table_data.append({
                    'Array': f"Array {i+1}",
                    'Accession': acc,
                    'Name': row.iloc[0]['name'],
                    'Organism': row.iloc[0]['organism_common'],
                    'Length': row.iloc[0]['sequence_length'],
                    'Function': function_display
                })
            else:
                table_data.append({
                    'Array': f"Array {i+1}",
                    'Accession': acc,
                    'Name': "Not found",
                    'Organism': "Not found",
                    'Length': "Not found",
                    'Function': "Not found"
                })
    return table_data

def create_kingdom_similarity_heatmap(df, accession_arrays, similarity_array, selected_kingdom):
    print(f"DEBUG: Creating sorted heatmap for kingdom: {selected_kingdom}")
    print(f"DEBUG: similarity_array shape: {similarity_array.shape if hasattr(similarity_array, 'shape') else 'not numpy array'}")
    print(f"DEBUG: accession_arrays length: {len(accession_arrays)}")
    
    # Create sorted arrays based on kingdom and similarity
    sorted_arrays = []
    sorted_similarities = []
    
    for i, accessions in enumerate(accession_arrays):
        accessions = [acc.decode('utf-8') if isinstance(acc, bytes) else acc for acc in accessions]
        if similarity_array is not None and len(similarity_array) > i:
            similarities = similarity_array[i].tolist() if hasattr(similarity_array[i], 'tolist') else similarity_array[i]
        else:
            similarities = [0.0] * len(accessions)
        
        # Get kingdom and similarity data for sorting
        protein_data = []
        for j, acc in enumerate(accessions):
            row = df[df['accession'] == acc]
            if not row.empty and pd.notna(row.iloc[0]['taxonomy_lineage']):
                taxonomy = row.iloc[0]['taxonomy_lineage']
                kingdom = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
            else:
                kingdom = 'Unknown'
            protein_data.append((acc, kingdom, similarities[j]))
        
        # Sort by kingdom first, then by similarity (descending) within each kingdom
        # Put selected kingdom first, then others
        def sort_key(item):
            acc, kingdom, sim = item
            if kingdom == selected_kingdom:
                return (0, -sim)  # Selected kingdom first, then by similarity (descending)
            else:
                return (1, kingdom, -sim)  # Other kingdoms, sorted by kingdom name, then similarity
        
        sorted_data = sorted(protein_data, key=sort_key)
        sorted_accessions = [acc for acc, kingdom, sim in sorted_data]
        sorted_sims = [sim for acc, kingdom, sim in sorted_data]
        
        sorted_arrays.append(sorted_accessions)
        sorted_similarities.append(sorted_sims)
    
    print(f"DEBUG: Sorted arrays: {sorted_arrays}")
    print(f"DEBUG: Sorted similarities: {sorted_similarities}")
    
    # Use the same heatmap creation logic as create_accession_heatmap
    max_len = max(len(arr) for arr in sorted_arrays)
    accession_labels = []
    kingdom_data = []
    hover_texts = []
    
    for i, accessions in enumerate(sorted_arrays):
        # Pad with empty strings to make all arrays the same length
        padded_accessions = accessions + [''] * (max_len - len(accessions))
        accession_labels.append(padded_accessions)
        
        # Get kingdom data for coloring
        row_kingdoms = []
        row_hover = []
        for j, acc in enumerate(accessions):
            row = df[df['accession'] == acc]
            if not row.empty and pd.notna(row.iloc[0]['taxonomy_lineage']):
                taxonomy = row.iloc[0]['taxonomy_lineage']
                kingdom = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
            else:
                kingdom = 'Unknown'
            row_kingdoms.append(kingdom)
            
            if not row.empty:
                function_text = row.iloc[0]['function']
                if pd.notna(function_text):
                    function_display = function_text[:50] + ('...' if len(function_text) > 50 else '')
                else:
                    function_display = 'No function data'
                taxonomy = row.iloc[0]['taxonomy_lineage'] if pd.notna(row.iloc[0]['taxonomy_lineage']) else 'Unknown'
                kingdom = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
                full_name = row.iloc[0]['full_name'] if pd.notna(row.iloc[0]['full_name']) else 'Not found'
                similarity = sorted_similarities[i][j]
                hover_text = f"""
                <b>Accession:</b> {acc}<br>
                <b>Full Name:</b> {full_name}<br>
                <b>Name:</b> {row.iloc[0]['name']}<br>
                <b>Organism:</b> {row.iloc[0]['organism_common']}<br>
                <b>Kingdom:</b> {kingdom}<br>
                <b>Similarity:</b> {similarity:.3f}<br>
                <b>Function:</b> {function_display}
                """
            else:
                hover_text = f"<b>Accession:</b> {acc}<br><b>Status:</b> Not found in database"
            row_hover.append(hover_text)
        
        # Pad kingdoms and hover with empty strings
        row_kingdoms.extend([''] * (max_len - len(accessions)))
        row_hover.extend([''] * (max_len - len(accessions)))
        kingdom_data.append(row_kingdoms)
        hover_texts.append(row_hover)
    
    # Transpose the data for the heatmap
    accession_labels_transposed = list(zip(*accession_labels))[::-1]
    hover_texts_transposed = list(zip(*hover_texts))[::-1]
    kingdom_data_transposed = list(zip(*kingdom_data))[::-1]
    
    # Create color mapping for kingdoms
    unique_kingdoms = sorted(set(k for row in kingdom_data_transposed for k in row if k and k != ''))
    kingdom_to_num = {kingdom: i for i, kingdom in enumerate(unique_kingdoms)}
    
    # Create numerical data based on kingdoms
    kingdom_numerical = []
    for row in kingdom_data_transposed:
        num_row = []
        for kingdom in row:
            if kingdom and kingdom != '':
                num_row.append(kingdom_to_num[kingdom])
            else:
                num_row.append(-1)  # -1 for empty cells
        kingdom_numerical.append(num_row)
    
    # Create color map for kingdoms
    kingdom_color_map = get_categorical_colormap(unique_kingdoms)
    
    # Custom colorscale for plotly
    colorscale = []
    for i, k in enumerate(unique_kingdoms):
        color = kingdom_color_map[k]
        if len(unique_kingdoms) == 1:
            colorscale = [[0, color], [1, color]]
        else:
            val = i / (len(unique_kingdoms) - 1)
            colorscale.append([val, color])
    if len(unique_kingdoms) > 1:
        last_color = kingdom_color_map[unique_kingdoms[-1]]
        colorscale.append([1.0, last_color])
    
    # Create the heatmap
    fig = go.Figure(data=go.Heatmap(
        z=kingdom_numerical,
        hovertext=hover_texts_transposed,
        hoverinfo='text',
        colorscale=colorscale,
        showscale=False,
        zmin=0,
        zmax=len(unique_kingdoms)-1
    ))
    
    # Add manual legend
    for i, k in enumerate(unique_kingdoms):
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(size=15, color=kingdom_color_map[k]),
            legendgroup=k,
            showlegend=True,
            name=k
        ))
    
    fig.update_layout(
        title=f"Sorted Heatmap: {selected_kingdom} first, then by similarity",
        xaxis_title="Array Index",
        yaxis_title="Position in Array (sorted by kingdom then similarity)",
        height=400,
        width=1200,
        legend_title="Kingdom"
    )
    
    print(f"DEBUG: Sorted heatmap created successfully")
    return fig

# --- Dash App ---
def run_dash_app(csv_path, accession_arrays, similarity_array=None):
    df = pd.read_csv(csv_path)
    kingdoms = get_available_kingdoms(df, accession_arrays)
    app = dash.Dash(__name__)
    app.layout = html.Div([
        html.H1("Protein Accession Visualization Dashboard (Dash)"),
        html.Label("Select Kingdom:"),
        dcc.Dropdown(
            id='kingdom-dropdown',
            options=[{'label': 'All Kingdoms', 'value': 'all'}] + [{'label': k, 'value': k} for k in kingdoms],
            value='all',
            clearable=False,
            style={'width': '300px'}
        ),
        dcc.Graph(id='heatmap-graph'),
        dcc.Graph(id='similarity-graph', style={'display': 'none'}),
        dcc.Graph(id='kingdom-sim-graph', style={'display': 'none'}),
        dash_table.DataTable(
            id='accession-table',
            columns=[{"name": i, "id": i} for i in ['Array', 'Accession', 'Name', 'Organism', 'Length', 'Function']],
            style_table={'overflowX': 'auto', 'width': '1200px'},
            style_cell={'textAlign': 'left', 'fontSize': 12},
            page_size=20
        )
    ])
    @app.callback(
        Output('heatmap-graph', 'figure'),
        Output('similarity-graph', 'figure'),
        Output('similarity-graph', 'style'),
        Output('kingdom-sim-graph', 'figure'),
        Output('kingdom-sim-graph', 'style'),
        Output('accession-table', 'data'),
        Input('kingdom-dropdown', 'value')
    )
    def update_dashboard(selected_kingdom):
        heatmap_fig = create_accession_heatmap(df, accession_arrays)
        sim_fig = go.Figure()
        sim_style = {'display': 'none'}
        kingdom_fig = go.Figure()
        kingdom_style = {'display': 'none'}
        if similarity_array is not None:
            sim_fig = create_similarity_heatmap(similarity_array)
            sim_style = {'display': 'block'}
            if selected_kingdom != 'all':
                kingdom_fig = create_kingdom_similarity_heatmap(df, accession_arrays, similarity_array, selected_kingdom)
                kingdom_style = {'display': 'block'}
        table_data = create_accession_table(df, accession_arrays)
        return heatmap_fig, sim_fig, sim_style, kingdom_fig, kingdom_style, table_data
    app.run(debug=True, host='0.0.0.0', port=8050)

# Example usage
if __name__ == "__main__":
    # accession_arrays = [
    #     ['P03756', 'Q21LK2', 'Q89EM9', 'Q03544', 'Q03546', 'P76515', 'A4IM80', 'B6JPK6', 'B5Z9N6', 'Q03547'],
    #     ['P43661', 'Q887Q8', 'Q8X5K6', 'P59791', 'P33128', 'Q88ND4', 'P70799', 'Q06062', 'P39632', 'P33409']
    # ]
    # similarity_array = [
    #     [1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.3, 0.5, 0.7, 0.9],
    #     [0.8, 1.0, 0.7, 0.5, 0.3, 0.2, 0.4, 0.6, 0.8, 0.7]
    # ]

    # load from numpy arrays
    accession_arrays = np.load("accession_arrays.npy")
    similarity_array = np.load("similarity_array.npy")
    print(f"DEBUG: accession_arrays shape: {accession_arrays.shape}")
    print(f"DEBUG: similarity_array shape: {similarity_array.shape}")
    print("Running Dash app...")
    csv_path = "/home/jovyan/workspace/data/unnotate/parsed_uniprot_swiss_data.csv"
    run_dash_app(csv_path, accession_arrays, similarity_array) 