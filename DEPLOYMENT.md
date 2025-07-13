# Deployment Guide

This guide explains how to deploy the Protein Accession Visualizer with automatic data loading from GitHub releases.

## Option 1: Streamlit Cloud (Recommended)

### Step 1: Push to GitHub
Ensure your repository is pushed to GitHub with the latest code.

### Step 2: Deploy to Streamlit Cloud
1. Go to [share.streamlit.io](https://share.streamlit.io)
2. Sign in with your GitHub account
3. Click "New app"
4. Select your repository: `khoa-yelo/unnotate`
5. Set the main file path to: `streamlit_app.py`
6. Click "Deploy"

### Step 3: Data Loading
The app will automatically:
1. Check for data files on startup
2. If missing, show a download interface
3. Download data from your GitHub release
4. Extract and load the data

## Option 2: Local Development

### Step 1: Clone and Setup
```bash
git clone https://github.com/khoa-yelo/unnotate.git
cd unnotate
pip install -r requirements.txt
```

### Step 2: Run the App
```bash
streamlit run streamlit_app.py
```

The app will automatically download data from GitHub releases if needed.

## Option 3: Manual Data Setup

If you prefer to manually manage data files:

1. **Download the data**:
   - Go to [GitHub Releases](https://github.com/khoa-yelo/unnotate/releases)
   - Download `protein_data_v1.0.zip`
   - Extract to your project directory

2. **Run the app**:
   ```bash
   streamlit run streamlit_app.py
   ```

## Data Files Required

The app needs these files in the project directory:
- `full_accession_arrays.npy` - Protein accession arrays
- `full_similarity_array.npy` - Similarity scores
- `parsed_uniprot_swiss_data.csv` - Protein metadata (680MB)

## Troubleshooting

### Data Download Issues
- Check your internet connection
- Verify the GitHub repository and release exist
- Try downloading manually from the releases page

### Memory Issues
- The CSV file is large (680MB)
- Ensure you have sufficient RAM (2GB+ recommended)
- Consider using a cloud deployment for better performance

### File Permission Issues
- Ensure the app has write permissions to extract downloaded files
- On Linux/Mac, you might need to run with appropriate permissions

## Environment Variables (Optional)

You can set these environment variables for customization:
- `GITHUB_REPO_OWNER` - Default repository owner
- `GITHUB_REPO_NAME` - Default repository name
- `GITHUB_RELEASE_TAG` - Default release tag

## Performance Tips

1. **First Run**: The initial data download may take a few minutes
2. **Subsequent Runs**: Data is cached and loads quickly
3. **Memory**: The app uses about 1-2GB RAM with full dataset
4. **Network**: Ensure stable internet for initial data download 