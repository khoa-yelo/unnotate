# Unnotate

## Overview
A tool designed to analyze hypothetical proteins. Unnotate generates protein embeddings using ESM-C and identifies their closest matches within a database of embeddings from reviewed Uniprot proteins. The findings can be interactively explored and analyzed through a Streamlit dashboard.

## Installation
```sh
git clone https://github.com/khoa-yelo/unnotate.git
```

### Using Conda/Mamba


#### GPU Option (Recommended if available - extremely fast)
```sh
mamba env create -f envs/unnotate_gpu.yaml
mamba activate unnotate_gpu
```

#### CPU Option
```sh
mamba env create -f envs/unnotate_cpu.yaml
mamba activate unnotate_cpu
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
  --fasta-file [your.faa] \ # protein sequences
  --database-dir [path/to/database/dir] \  # dir obtained from download_db
  --output-dir [path/to/output/dir] \
  --k 20 # number of nearest neighbors
  --prefix [mysterious_virus]
  --cpu # specify if using cpu env, remove if using gpu 
```

#### Outputs
- `PREFIX_uniprot.csv`: uniprot information of the accession hits (nearest neighbors from query fasta)
- `PREFIX_results.npz`: numpy zip file containing 
    - accession: Uniprot accesssion nearest neighbors of queried proteins (Nxk)
    - cosine_similarity: ESM-C embeddings cosine similarity between nearest neighbors and queried proteins (Nxk)
    - sequence_identity: pairwise sequence identity between nearest neighbors and queried proteins (Nxk)
- `PREFIX_streamlit.zip`: zip file containing `PREFIX_uniprot.csv`, `PREFIX_results.npz` for visualization
- statistics (coming soon!)

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