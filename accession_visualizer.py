"""
Interactive Accession Visualizer
"""

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
from typing import List, Union, Optional
import warnings
warnings.filterwarnings('ignore')

# Helper for categorical color mapping
from matplotlib import cm, colors as mcolors

def get_categorical_colormap(categories):
    # Use tab20 or similar for up to 20 categories
    base_colors = list(mcolors.TABLEAU_COLORS.values()) + list(mcolors.CSS4_COLORS.values())
    color_map = {}
    for i, cat in enumerate(categories):
        color_map[cat] = base_colors[i % len(base_colors)]
    return color_map

class AccessionVisualizer:
    """
    A class for creating interactive visualizations of protein accession arrays.
    """
    
    def __init__(self, csv_path: str):
        """
        Initialize the visualizer with the CSV data.
        
        Args:
            csv_path (str): Path to the parsed_uniprot_swiss_data.csv file.
        """
        self.csv_path = csv_path
        self.df = pd.read_csv(csv_path)
        print(f"Loaded {len(self.df)} protein records")
        
    def create_accession_grid(self, accession_arrays: List[List[str]], 
                            title: str = "Protein Accession Visualization",
                            colors: List[str] = None) -> go.Figure:
        """
        Create an interactive grid visualization of accession arrays.
        
        Args:
            accession_arrays (List[List[str]]): List of accession arrays to visualize.
            title (str): Title for the visualization.
            colors (List[str]): List of colors for each array.
            
        Returns:
            go.Figure: Interactive plotly figure.
        """
        if colors is None:
            colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                     '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        
        # Create subplots
        n_arrays = len(accession_arrays)
        
        # Adjust spacing for large numbers of arrays
        if n_arrays > 10:
            vertical_spacing = 0.01
            height_per_row = 100
        else:
            vertical_spacing = 0.1
            height_per_row = 200
        
        fig = make_subplots(
            rows=n_arrays, cols=1,
            subplot_titles=[f"Array {i+1}" for i in range(n_arrays)],
            vertical_spacing=vertical_spacing
        )
        
        for i, accessions in enumerate(accession_arrays):
            # Convert bytes to strings if needed
            accessions = [acc.decode('utf-8') if isinstance(acc, bytes) else acc for acc in accessions]
            
            # Get data for these accessions
            accession_data = []
            for acc in accessions:
                row = self.df[self.df['accession'] == acc]
                if not row.empty:
                    accession_data.append({
                        'accession': acc,
                        'name': row.iloc[0]['name'],
                        'full_name': row.iloc[0]['full_name'],
                        'organism': row.iloc[0]['organism'],
                        'organism_common': row.iloc[0]['organism_common'],
                        'sequence_length': row.iloc[0]['sequence_length'],
                        'function': row.iloc[0]['function'] if pd.notna(row.iloc[0]['function']) else 'No function data',
                        'protein_families': row.iloc[0]['protein_families'] if pd.notna(row.iloc[0]['protein_families']) else 'No family data'
                    })
                else:
                    accession_data.append({
                        'accession': acc,
                        'name': 'Not found',
                        'full_name': 'Not found',
                        'organism': 'Not found',
                        'organism_common': 'Not found',
                        'sequence_length': 'Not found',
                        'function': 'Not found',
                        'protein_families': 'Not found'
                    })
            
            # Create hover text
            hover_texts = []
            for data in accession_data:
                hover_text = f"""
                <b>Accession:</b> {data['accession']}<br>
                <b>Name:</b> {data['name']}<br>
                <b>Full Name:</b> {data['full_name']}<br>
                <b>Organism:</b> {data['organism']}<br>
                <b>Common Name:</b> {data['organism_common']}<br>
                <b>Sequence Length:</b> {data['sequence_length']}<br>
                <b>Function:</b> {data['function'][:100]}{'...' if len(data['function']) > 100 else ''}<br>
                <b>Protein Families:</b> {data['protein_families'][:100]}{'...' if len(data['protein_families']) > 100 else ''}
                """
                hover_texts.append(hover_text)
            
            # Create the grid
            x_positions = list(range(len(accessions)))
            y_positions = [i] * len(accessions)
            
            fig.add_trace(
                go.Scatter(
                    x=x_positions,
                    y=y_positions,
                    mode='markers+text',
                    marker=dict(
                        size=20,
                        color=colors[i % len(colors)],
                        symbol='square',
                        line=dict(width=2, color='black')
                    ),
                    text=accessions,
                    textposition="middle center",
                    textfont=dict(size=10, color='white'),
                    hovertext=hover_texts,
                    hoverinfo='text',
                    name=f'Array {i+1}',
                    showlegend=False
                ),
                row=i+1, col=1
            )
            
            # Update x-axis for each subplot
            fig.update_xaxes(
                range=[-0.5, len(accessions) - 0.5],
                showticklabels=False,
                row=i+1, col=1
            )
            
            # Update y-axis for each subplot
            fig.update_yaxes(
                range=[i-0.3, i+0.3],
                showticklabels=False,
                row=i+1, col=1
            )
        
        # Update layout
        fig.update_layout(
            title=title,
            height=height_per_row * n_arrays,
            width=1200,
            showlegend=False,
            plot_bgcolor='white'
        )
        
        return fig
    
    def create_accession_heatmap(self, accession_arrays: List[List[str]], 
                               title: str = "Protein Accession Heatmap") -> go.Figure:
        """
        Create a heatmap visualization of accession arrays.
        
        Args:
            accession_arrays (List[List[str]]): List of accession arrays to visualize.
            title (str): Title for the visualization.
            
        Returns:
            go.Figure: Interactive plotly figure.
        """
        # Prepare data for heatmap
        max_len = max(len(arr) for arr in accession_arrays)
        heatmap_data = []
        accession_labels = []
        kingdom_data = []
        full_name_data = []
        
        for i, accessions in enumerate(accession_arrays):
            # Convert bytes to strings if needed
            accessions = [acc.decode('utf-8') if isinstance(acc, bytes) else acc for acc in accessions]
            
            # Pad with empty strings to make all arrays the same length
            padded_accessions = accessions + [''] * (max_len - len(accessions))
            accession_labels.append(padded_accessions)
            
            # Get kingdom data for coloring
            row_kingdoms = []
            row_full_names = []
            for acc in accessions:
                row = self.df[self.df['accession'] == acc]
                if not row.empty and pd.notna(row.iloc[0]['taxonomy_lineage']):
                    taxonomy = row.iloc[0]['taxonomy_lineage']
                    kingdom = taxonomy.split(';')[0].strip() if ';' in taxonomy else taxonomy.strip()
                else:
                    kingdom = 'Unknown'
                row_kingdoms.append(kingdom)
                if not row.empty and pd.notna(row.iloc[0]['full_name']):
                    full_name = row.iloc[0]['full_name']
                else:
                    full_name = 'Not found'
                row_full_names.append(full_name)
            
            # Pad kingdoms with empty strings
            row_kingdoms.extend([''] * (max_len - len(accessions)))
            row_full_names.extend([''] * (max_len - len(accessions)))
            kingdom_data.append(row_kingdoms)
            full_name_data.append(row_full_names)
            
            # Create numerical data for heatmap (using array index as value)
            row_data = [i] * len(padded_accessions)
            heatmap_data.append(row_data)
        
        # Create hover text
        hover_texts = []
        for i, accessions in enumerate(accession_arrays):
            accessions = [acc.decode('utf-8') if isinstance(acc, bytes) else acc for acc in accessions]
            row_hover = []
            
            for j, acc in enumerate(accessions):
                row = self.df[self.df['accession'] == acc]
                if not row.empty:
                    function_text = row.iloc[0]['function']
                    if pd.notna(function_text):
                        function_display = function_text[:50] + ('...' if len(function_text) > 50 else '')
                    else:
                        function_display = 'No function data'
                    
                    # Get taxonomy info
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
            
            # Pad with empty hover text
            row_hover.extend([''] * (max_len - len(accessions)))
            hover_texts.append(row_hover)
        
        # Transpose the data for the heatmap
        accession_labels_transposed = list(zip(*accession_labels))
        hover_texts_transposed = list(zip(*hover_texts))
        kingdom_data_transposed = list(zip(*kingdom_data))
        full_name_data_transposed = list(zip(*full_name_data))
        
        # Flip vertically (reverse the order so index 0 is at top)
        accession_labels_transposed = accession_labels_transposed[::-1]
        hover_texts_transposed = hover_texts_transposed[::-1]
        kingdom_data_transposed = kingdom_data_transposed[::-1]
        full_name_data_transposed = full_name_data_transposed[::-1]
        
        # Create color mapping for kingdoms
        unique_kingdoms = sorted(set(k for row in kingdom_data_transposed for k in row if k and k != ''))
        
        # Create numerical mapping for kingdoms
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
        
        # Create color mapping for full names
        unique_full_names = sorted(set(fn for row in full_name_data_transposed for fn in row if fn and fn != ''))
        full_name_to_num = {fn: i for i, fn in enumerate(unique_full_names)}
        
        # Create numerical data based on full names
        full_name_numerical = []
        for row in full_name_data_transposed:
            num_row = []
            for full_name in row:
                if full_name and full_name != '':
                    num_row.append(full_name_to_num[full_name])
                else:
                    num_row.append(-1) # -1 for empty cells
            full_name_numerical.append(num_row)
        
        # Create color map for all kingdoms at once
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
        # Add the last color at position 1.0
        if len(unique_kingdoms) > 1:
            last_color = kingdom_color_map[unique_kingdoms[-1]]
            colorscale.append([1.0, last_color])
        
        # Plotly Heatmap
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
            title=title,
            xaxis_title="Array Index",
            yaxis_title="Position in Array",
            height=400,
            width=1200,
            legend_title="Kingdom"
        )
        
        return fig
    
    def create_similarity_heatmap(self, similarity_array: List[List[float]], title: str = "Similarity Heatmap") -> go.Figure:
        """
        Create a heatmap visualization of accession arrays.
        
        Args:
            accession_arrays (List[List[str]]): List of accession arrays to visualize.
            title (str): Title for the visualization.
            
        Returns:
            go.Figure: Interactive plotly figure.
        """
        # Transpose and flip vertically
        z = np.array(similarity_array).T
        z = z[::-1]  # Flip vertically so index 0 is at top
        
        fig = go.Figure(data=go.Heatmap(
            z=z,
            colorscale='Blues',
            colorbar=dict(title="Similarity")
        ))
        fig.update_layout(
            title=title,
            xaxis_title="Array Index",
            yaxis_title="Position in Array",
            height=400,
            width=1200
        )
        return fig

    def create_accession_table(self, accession_arrays: List[List[str]], 
                             title: str = "Protein Accession Information") -> go.Figure:
        """
        Create an interactive table visualization of accession information.
        
        Args:
            accession_arrays (List[List[str]]): List of accession arrays to visualize.
            title (str): Title for the visualization.
            
        Returns:
            go.Figure: Interactive plotly figure.
        """
        # Prepare table data
        table_data = []
        
        for i, accessions in enumerate(accession_arrays):
            accessions = [acc.decode('utf-8') if isinstance(acc, bytes) else acc for acc in accessions]
            
            for j, acc in enumerate(accessions):
                row = self.df[self.df['accession'] == acc]
                if not row.empty:
                    function_text = row.iloc[0]['function']
                    if pd.notna(function_text):
                        function_display = function_text[:50] + '...' if len(function_text) > 50 else function_text
                    else:
                        function_display = 'No function data'
                    
                    table_data.append([
                        f"Array {i+1}",
                        acc,
                        row.iloc[0]['name'],
                        row.iloc[0]['organism_common'],
                        row.iloc[0]['sequence_length'],
                        function_display
                    ])
                else:
                    table_data.append([
                        f"Array {i+1}",
                        acc,
                        "Not found",
                        "Not found",
                        "Not found",
                        "Not found"
                    ])
        
        fig = go.Figure(data=[go.Table(
            header=dict(
                values=['Array', 'Accession', 'Name', 'Organism', 'Length', 'Function'],
                fill_color='paleturquoise',
                align='left',
                font=dict(size=12)
            ),
            cells=dict(
                values=list(zip(*table_data)),
                fill_color='lavender',
                align='left',
                font=dict(size=10)
            )
        )])
        
        fig.update_layout(
            title=title,
            height=600,
            width=1200
        )
        
        return fig


def visualize_accessions(csv_path: str, accession_arrays: List[List[str]], output_file: Optional[str] = None, show_plot: bool = True, similarity_array: Optional[List[List[float]]] = None):
    """
    Convenience function to create and display accession visualizations.
    
    Args:
        csv_path (str): Path to the CSV file.
        accession_arrays (List[List[str]]): List of accession arrays to visualize.
        output_file (str, optional): Path to save the HTML file.
        show_plot (bool): Whether to display the plot.
    """
    visualizer = AccessionVisualizer(csv_path)
    
    # Create visualizations (removed grid)
    heatmap_fig = visualizer.create_accession_heatmap(accession_arrays)
    figs = [heatmap_fig]
    subplot_titles = ["Kingdom Heatmap"]
    specs = [[{"type": "heatmap"}]]
    rows = 1
    if similarity_array is not None:
        sim_fig = visualizer.create_similarity_heatmap(similarity_array)
        figs.append(sim_fig)
        subplot_titles.append("Similarity Heatmap")
        specs.append([{"type": "heatmap"}])
        rows += 1
    table_fig = visualizer.create_accession_table(accession_arrays)
    figs.append(table_fig)
    subplot_titles.append("Table View")
    specs.append([{"type": "table"}])
    rows += 1
    combined_fig = make_subplots(
        rows=rows, cols=1,
        subplot_titles=subplot_titles,
        vertical_spacing=0.1,
        specs=specs
    )
    
    for i, fig in enumerate(figs):
        for trace in fig.data:
            combined_fig.add_trace(trace, row=i+1, col=1)
    
    # Update layout
    n_arrays = len(accession_arrays)
    if n_arrays > 10:
        dashboard_height = 600
    else:
        dashboard_height = 800
        
    combined_fig.update_layout(
        title="Protein Accession Visualization Dashboard",
        height=dashboard_height,
        width=1200,
        showlegend=True
    )
    
    if output_file:
        combined_fig.write_html(output_file)
        print(f"Visualization saved to {output_file}")
    
    if show_plot:
        try:
            combined_fig.show()
        except Exception as e:
            print(f"Could not display plot: {e}")
            print("Plot saved to file instead.")
    
    return combined_fig


# Example usage
if __name__ == "__main__":
    # Example accession arrays
    accession_arrays = [
        ['P03756', 'Q21LK2', 'Q89EM9', 'Q03544', 'Q03546', 'P76515', 'A4IM80', 'B6JPK6', 'B5Z9N6', 'Q03547'],
        ['P43661', 'Q887Q8', 'Q8X5K6', 'P59791', 'P33128', 'Q88ND4', 'P70799', 'Q06062', 'P39632', 'P33409']
    ]
    
    # Example similarity array
    similarity_array = [
        [1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.3, 0.5, 0.7, 0.9],
        [0.8, 1.0, 0.7, 0.5, 0.3, 0.2, 0.4, 0.6, 0.8, 0.7]
    ]
    
    csv_path = "/home/jovyan/workspace/data/unnotate/parsed_uniprot_swiss_data.csv"
    
    # Create visualization
    fig = visualize_accessions(csv_path, accession_arrays, similarity_array=similarity_array, 
                             output_file="accession_visualization.html", show_plot=False)
