# Unnotate

## Overview
Unnotate is a tool to annotate proteins that traditional annotation pipelines miss. It works by creating protein embeddings with the ESM-C model and comparing them against a comprehensive database of embeddings from Uniprot's reviewed proteins. The results can be explored through an interactive Streamlit dashboard for easy visualization and analysis.

## Installation

### Using Conda/Mamba

#### GPU Option (extremely fast)
```sh
git clone https://github.com/khoa-yelo/unnotate.git
mamba env create -f envs/unnotate_gpu.yaml
mamba activate unnotate_gpu
```
## Usage

### CLI
After installation, use the `unnotate` bash script:

#### Download the database
```sh
./unnotate download_db --dest ./database
```

#### Annotate protein sequences
```sh
./unnotate unnot \
  --fasta-file [your.faa] \  # protein sequences
  --database-dir [path/to/database_folder] \  # folder containing .h5 and .csv
  --output-dir ./results \
  --k 20 # number of nearest neighbors
  --prefix [mysterious_virus]
```

## Visualize Results with Streamlit

After running annotation, you can explore your results interactively using the Unnotate Streamlit dashboard.

### Option 1: Run Streamlit Locally

1. **Start the dashboard:**
   ```sh
   streamlit run streamlit_app.py
   ```
2. **Open your browser** to the address shown in the terminal (usually http://localhost:8501).
3. **Upload your results:**
   - Click the "Upload Data from ZIP File" expander.
   - Upload the output zip file from the annotation step (e.g., `mysterious_virus_streamlit.zip`).

### Option 2: Use the Public Streamlit Cloud

- Go to [https://unnotate.streamlit.app](https://unnotate.streamlit.app)
- Upload your output zip file as above.

The Streamlit dashboard provides access to three pre-loaded example datasets in addition to any custom data you upload: Virus, Bacteria, and Putative Phage Plasmid datasets.