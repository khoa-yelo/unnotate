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

class FastAccessionVisualizer:
    def __init__(self, csv_path, accession_arrays, similarity_array=None):
        """Initialize with pre-processed data for fast lookups"""
        print("Loading and pre-processing data...")
        self.df = pd.read_csv(csv_path)
        
        # Convert accession arrays to list format and handle bytes
        self.accession_arrays = []
        for arr in accession_arrays:
            if hasattr(arr, 'tolist'):
                arr = arr.tolist()
            arr = [acc.decode('utf-8') if isinstance(acc, bytes) else acc for acc in arr]
            self.accession_arrays.append(arr)
        
        self.similarity_array = similarity_array
        
        # Create fast lookup dictionary
        self.accession_data = {}
        for _, row in self.df.iterrows():
            acc = row['accession']
            if pd.notna(row['taxonomy_lineage']):
                taxonomy = row['taxonomy_lineage']
                kingdom = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
            else:
                kingdom = 'Unknown'
            
            self.accession_data[acc] = {
                'name': row['name'],
                'full_name': row['full_name'] if pd.notna(row['full_name']) else 'Not found',
                'organism': row['organism'],
                'organism_common': row['organism_common'],
                'sequence_length': row['sequence_length'],
                'function': row['function'] if pd.notna(row['function']) else 'No function data',
                'kingdom': kingdom
            }
        
        # Pre-calculate all kingdoms
        self.kingdoms = sorted(set(data['kingdom'] for data in self.accession_data.values()))
        
        # Pre-calculate all figures
        print("Pre-calculating figures...")
        self.heatmap_fig = self._create_accession_heatmap()
        self.sim_fig = self._create_similarity_heatmap() if similarity_array is not None else go.Figure()
        self.kingdom_figs = {}
        for kingdom in self.kingdoms:
            self.kingdom_figs[kingdom] = self._create_kingdom_similarity_heatmap(kingdom)
        
        # Pre-calculate table data
        self.table_data = self._create_accession_table()
        
        print("Data pre-processing complete!")

    def _create_accession_heatmap(self):
        """Create kingdom heatmap with pre-processed data"""
        max_len = max(len(arr) for arr in self.accession_arrays)
        kingdom_data = []
        hover_texts = []
        
        for accessions in self.accession_arrays:
            row_kingdoms = []
            row_hover = []
            
            for acc in accessions:
                if acc in self.accession_data:
                    data = self.accession_data[acc]
                    row_kingdoms.append(data['kingdom'])
                    
                    function_display = data['function'][:50] + ('...' if len(data['function']) > 50 else '')
                    hover_text = f"""
                    <b>Accession:</b> {acc}<br>
                    <b>Full Name:</b> {data['full_name']}<br>
                    <b>Name:</b> {data['name']}<br>
                    <b>Organism:</b> {data['organism_common']}<br>
                    <b>Kingdom:</b> {data['kingdom']}<br>
                    <b>Function:</b> {function_display}
                    """
                else:
                    row_kingdoms.append('Unknown')
                    hover_text = f"<b>Accession:</b> {acc}<br><b>Status:</b> Not found in database"
                row_hover.append(hover_text)
            
            # Pad with empty strings
            row_kingdoms.extend([''] * (max_len - len(accessions)))
            row_hover.extend([''] * (max_len - len(accessions)))
            kingdom_data.append(row_kingdoms)
            hover_texts.append(row_hover)
        
        # Transpose and flip
        kingdom_data_transposed = list(zip(*kingdom_data))[::-1]
        hover_texts_transposed = list(zip(*hover_texts))[::-1]
        
        # Create numerical mapping
        kingdom_to_num = {kingdom: i for i, kingdom in enumerate(self.kingdoms)}
        kingdom_numerical = []
        for row in kingdom_data_transposed:
            num_row = []
            for kingdom in row:
                if kingdom and kingdom != '':
                    num_row.append(kingdom_to_num[kingdom])
                else:
                    num_row.append(-1)
            kingdom_numerical.append(num_row)
        
        # Create color map
        kingdom_color_map = get_categorical_colormap(self.kingdoms)
        colorscale = []
        for i, k in enumerate(self.kingdoms):
            color = kingdom_color_map[k]
            if len(self.kingdoms) == 1:
                colorscale = [[0, color], [1, color]]
            else:
                val = i / (len(self.kingdoms) - 1)
                colorscale.append([val, color])
        if len(self.kingdoms) > 1:
            last_color = kingdom_color_map[self.kingdoms[-1]]
            colorscale.append([1.0, last_color])
        
        # Create figure
        fig = go.Figure(data=go.Heatmap(
            z=kingdom_numerical,
            hovertext=hover_texts_transposed,
            hoverinfo='text',
            colorscale=colorscale,
            showscale=False,
            zmin=0,
            zmax=len(self.kingdoms)-1
        ))
        
        # Add legend
        for k in self.kingdoms:
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

    def _create_similarity_heatmap(self):
        """Create similarity heatmap"""
        z = np.array(self.similarity_array).T[::-1]
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

    def _create_kingdom_similarity_heatmap(self, selected_kingdom):
        """Create kingdom-specific sorted heatmap"""
        # Create sorted arrays based on kingdom and similarity
        sorted_arrays = []
        sorted_similarities = []
        
        for i, accessions in enumerate(self.accession_arrays):
            similarities = self.similarity_array[i].tolist() if self.similarity_array is not None else [0.0] * len(accessions)
            
            # Get protein data for sorting
            protein_data = []
            for j, acc in enumerate(accessions):
                if acc in self.accession_data:
                    kingdom = self.accession_data[acc]['kingdom']
                else:
                    kingdom = 'Unknown'
                protein_data.append((acc, kingdom, similarities[j]))
            
            # Sort by kingdom first, then by similarity
            def sort_key(item):
                acc, kingdom, sim = item
                if kingdom == selected_kingdom:
                    return (0, -sim)
                else:
                    return (1, kingdom, -sim)
            
            sorted_data = sorted(protein_data, key=sort_key)
            sorted_accessions = [acc for acc, kingdom, sim in sorted_data]
            sorted_sims = [sim for acc, kingdom, sim in sorted_data]
            
            sorted_arrays.append(sorted_accessions)
            sorted_similarities.append(sorted_sims)
        
        # Create heatmap using sorted data
        max_len = max(len(arr) for arr in sorted_arrays)
        kingdom_data = []
        hover_texts = []
        
        for i, accessions in enumerate(sorted_arrays):
            row_kingdoms = []
            row_hover = []
            
            for j, acc in enumerate(accessions):
                if acc in self.accession_data:
                    data = self.accession_data[acc]
                    row_kingdoms.append(data['kingdom'])
                    
                    function_display = data['function'][:50] + ('...' if len(data['function']) > 50 else '')
                    similarity = sorted_similarities[i][j]
                    hover_text = f"""
                    <b>Accession:</b> {acc}<br>
                    <b>Full Name:</b> {data['full_name']}<br>
                    <b>Name:</b> {data['name']}<br>
                    <b>Organism:</b> {data['organism_common']}<br>
                    <b>Kingdom:</b> {data['kingdom']}<br>
                    <b>Similarity:</b> {similarity:.3f}<br>
                    <b>Function:</b> {function_display}
                    """
                else:
                    row_kingdoms.append('Unknown')
                    hover_text = f"<b>Accession:</b> {acc}<br><b>Status:</b> Not found in database"
                row_hover.append(hover_text)
            
            # Pad with empty strings
            row_kingdoms.extend([''] * (max_len - len(accessions)))
            row_hover.extend([''] * (max_len - len(accessions)))
            kingdom_data.append(row_kingdoms)
            hover_texts.append(row_hover)
        
        # Transpose and flip
        kingdom_data_transposed = list(zip(*kingdom_data))[::-1]
        hover_texts_transposed = list(zip(*hover_texts))[::-1]
        
        # Create numerical mapping
        kingdom_to_num = {kingdom: i for i, kingdom in enumerate(self.kingdoms)}
        kingdom_numerical = []
        for row in kingdom_data_transposed:
            num_row = []
            for kingdom in row:
                if kingdom and kingdom != '':
                    num_row.append(kingdom_to_num[kingdom])
                else:
                    num_row.append(-1)
            kingdom_numerical.append(num_row)
        
        # Create color map
        kingdom_color_map = get_categorical_colormap(self.kingdoms)
        colorscale = []
        for i, k in enumerate(self.kingdoms):
            color = kingdom_color_map[k]
            if len(self.kingdoms) == 1:
                colorscale = [[0, color], [1, color]]
            else:
                val = i / (len(self.kingdoms) - 1)
                colorscale.append([val, color])
        if len(self.kingdoms) > 1:
            last_color = kingdom_color_map[self.kingdoms[-1]]
            colorscale.append([1.0, last_color])
        
        # Create figure
        fig = go.Figure(data=go.Heatmap(
            z=kingdom_numerical,
            hovertext=hover_texts_transposed,
            hoverinfo='text',
            colorscale=colorscale,
            showscale=False,
            zmin=0,
            zmax=len(self.kingdoms)-1
        ))
        
        # Add legend
        for k in self.kingdoms:
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
        return fig

    def _create_accession_table(self):
        """Create table data with pre-processed data"""
        table_data = []
        for i, accessions in enumerate(self.accession_arrays):
            for acc in accessions:
                if acc in self.accession_data:
                    data = self.accession_data[acc]
                    function_display = data['function'][:50] + '...' if len(data['function']) > 50 else data['function']
                    table_data.append({
                        'Array': f"Array {i+1}",
                        'Accession': acc,
                        'Name': data['name'],
                        'Organism': data['organism_common'],
                        'Length': data['sequence_length'],
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

def run_fast_dash_app(csv_path, accession_arrays, similarity_array=None):
    """Run optimized Dash app with pre-processed data"""
    # Initialize the fast visualizer
    visualizer = FastAccessionVisualizer(csv_path, accession_arrays, similarity_array)
    
    # Create Dash app
    app = dash.Dash(__name__)
    app.layout = html.Div([
        html.H1("Protein Accession Visualization Dashboard (Fast)"),
        html.Label("Select Kingdom:"),
        dcc.Dropdown(
            id='kingdom-dropdown',
            options=[{'label': 'All Kingdoms', 'value': 'all'}] + [{'label': k, 'value': k} for k in visualizer.kingdoms],
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
        # Return pre-calculated figures - no processing needed!
        sim_style = {'display': 'block'} if visualizer.similarity_array is not None else {'display': 'none'}
        kingdom_style = {'display': 'block'} if selected_kingdom != 'all' else {'display': 'none'}
        kingdom_fig = visualizer.kingdom_figs.get(selected_kingdom, go.Figure()) if selected_kingdom != 'all' else go.Figure()
        
        return (
            visualizer.heatmap_fig,
            visualizer.sim_fig,
            sim_style,
            kingdom_fig,
            kingdom_style,
            visualizer.table_data
        )
    
    app.run(debug=True, host='0.0.0.0', port=8050)

# Example usage
if __name__ == "__main__":
    # Test with example data to demonstrate speed improvement
    accession_arrays = [
        ['P03756', 'Q21LK2', 'Q89EM9', 'Q03544', 'Q03546', 'P76515', 'A4IM80', 'B6JPK6', 'B5Z9N6', 'Q03547'],
        ['P43661', 'Q887Q8', 'Q8X5K6', 'P59791', 'P33128', 'Q88ND4', 'P70799', 'Q06062', 'P39632', 'P33409']
    ]
    similarity_array = np.array([
        [1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.3, 0.5, 0.7, 0.9],
        [0.8, 1.0, 0.7, 0.5, 0.3, 0.2, 0.4, 0.6, 0.8, 0.7]
    ])
    # load from numpy arrays
    accession_arrays = np.load("accession_arrays.npy")
    similarity_array = np.load("similarity_array.npy")
    print(f"DEBUG: accession_arrays length: {len(accession_arrays)}")
    print(f"DEBUG: similarity_array shape: {similarity_array.shape}")
    print("Running Fast Dash app...")
    csv_path = "/home/jovyan/workspace/data/unnotate/parsed_uniprot_swiss_data.csv"
    run_fast_dash_app(csv_path, accession_arrays, similarity_array) 