# Protein Accession Visualizer

A Streamlit application for visualizing protein accession data with interactive heatmaps and detailed protein information tables.

## Features

- **Interactive Heatmaps**: Domain-based and similarity heatmaps with hover information
- **Domain-based Sorting**: Sort proteins by taxonomic domain and similarity scores
- **Protein Information Table**: Detailed table with clickable UniProt links
- **Flexible Data Loading**: Support for multiple data sources and file formats
- **Performance Optimized**: Efficient batch processing and vectorized operations

## Quick Start

### Option 1: Download from GitHub Release (Recommended)

1. **Download the data files**:
   - Go to the [GitHub Releases page](https://github.com/khoa-yelo/unnotate/releases)
   - Download the latest `protein_data_v1.0.zip` file
   - Extract the zip file to your project directory

2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Run the app**:
   ```bash
   streamlit run streamlit_app.py
   ```

### Option 2: Upload Your Own Data

1. **Prepare your data files**:
   - `accession_arrays.npy` - Protein accession arrays
   - `similarity_array.npy` - Similarity scores
   - `parsed_uniprot_swiss_data.csv` - Protein metadata

2. **Run the data downloader**:
   ```bash
   streamlit run data_downloader.py
   ```

3. **Upload your files** using the web interface

4. **Run the main app**:
   ```bash
   streamlit run streamlit_app.py
   ```

## Data Files

The application requires the following data files:

- **`full_accession_arrays.npy`**: Array of protein accession numbers organized by CDS
- **`full_similarity_array.npy`**: Similarity scores for protein comparisons
- **`parsed_uniprot_swiss_data.csv`**: Protein metadata from UniProt including:
  - Accession numbers
  - Protein names and full names
  - Organism information
  - Taxonomic classification
  - Function descriptions
  - Sequence lengths

## Usage

1. **Select a Domain**: Use the sidebar to highlight proteins from specific taxonomic domains
2. **View Heatmaps**: 
   - Domain Heatmap: Shows proteins colored by taxonomic domain
   - Similarity Heatmap: Shows similarity scores between proteins
3. **Explore Details**: Click on accession numbers in the table to open UniProt pages
4. **Sort Data**: When a domain is selected, both heatmaps and table are sorted accordingly

## Deployment

### Streamlit Cloud

1. **Push to GitHub**: Ensure your repository is on GitHub
2. **Deploy**: Connect your repository to [Streamlit Cloud](https://streamlit.io/cloud)
3. **Configure**: Set the main file path to `streamlit_app.py`

### Local Deployment

```bash
# Install dependencies
pip install -r requirements.txt

# Run the app
streamlit run streamlit_app.py
```

## File Structure

```
unnotate/
├── streamlit_app.py          # Main Streamlit application
├── data_downloader.py        # Data file management utility
├── create_release.py         # GitHub release creation script
├── requirements.txt          # Python dependencies
├── .streamlit/
│   └── config.toml          # Streamlit configuration
├── README.md                # This file
└── data files (npy, csv)    # Your protein data files
```

## Creating New Releases

To create a new data release:

```bash
python create_release.py
```

This will:
1. Create a zip file with all data files
2. Create a git tag
3. Push the tag to GitHub
4. Provide instructions for creating the GitHub release

## Support

For issues or questions:
1. Check the data files are in the correct format
2. Ensure all dependencies are installed
3. Try the data downloader utility for file management
