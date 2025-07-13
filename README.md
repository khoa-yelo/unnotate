# Protein Accession Visualizer

A Streamlit web application for visualizing protein accession data with interactive heatmaps and domain-based sorting.

## Features

- **Interactive Domain Heatmap**: Visualize protein domains with color-coded kingdoms
- **Similarity Heatmap**: View similarity scores between proteins
- **Domain-based Sorting**: Sort proteins by selected domain and similarity
- **Interactive Table**: Browse protein information with clickable UniProt links
- **Performance Optimized**: Uses vectorized pandas operations for fast data processing

## Installation

1. Clone the repository:
```bash
git clone <your-repo-url>
cd unnotate
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Local Development

Run the Streamlit app locally:
```bash
streamlit run streamlit_app.py
```

### Streamlit Cloud Deployment

1. Push your code to GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect your GitHub repository
4. Set the main file path to: `unnotate/streamlit_app.py`
5. Deploy!

## Data Requirements

The app expects:
- A CSV file with protein data (columns: accession, name, full_name, organism_common, taxonomy_lineage, sequence_length, function)
- NumPy arrays for accession arrays and similarity data

## Configuration

- Modify data file paths in `streamlit_app.py`
- Adjust visualization settings in the app
- Customize color schemes and layout as needed

## License

[Your License Here]
