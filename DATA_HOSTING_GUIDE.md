# Data Hosting Guide for Protein Accession Visualizer

This guide covers different options for hosting your data files for the Streamlit app deployment.

## 📊 Your Data Files

Based on your current setup, you have these data files:
- `full_accession_arrays.npy` (~63KB)
- `full_similarity_array.npy` (~25KB) 
- `parsed_uniprot_swiss_data.csv` (protein database)
- `embeddings.h5` (~9.9GB - very large!)

## 🚀 Hosting Options

### Option 1: GitHub Releases (Recommended for < 100MB files)

**Best for**: Small to medium files (< 100MB)

**Steps**:
1. Go to your GitHub repository: `https://github.com/khoa-yelo/unnotate`
2. Click "Releases" → "Create a new release"
3. Tag: `v1.0.0`
4. Title: `Initial Data Release`
5. Upload your files:
   - `full_accession_arrays.npy`
   - `full_similarity_array.npy`
   - `parsed_uniprot_swiss_data.csv`
6. Publish release

**Usage in app**:
```python
# URLs will be:
# https://github.com/khoa-yelo/unnotate/releases/download/v1.0.0/full_accession_arrays.npy
# https://github.com/khoa-yelo/unnotate/releases/download/v1.0.0/full_similarity_array.npy
# https://github.com/khoa-yelo/unnotate/releases/download/v1.0.0/parsed_uniprot_swiss_data.csv
```

### Option 2: Google Drive (Recommended for large files)

**Best for**: Large files like `embeddings.h5`

**Steps**:
1. Upload files to Google Drive
2. Right-click → "Get shareable link"
3. Convert to direct download link:
   - Replace `https://drive.google.com/file/d/FILE_ID/view?usp=sharing`
   - With `https://drive.google.com/uc?id=FILE_ID&export=download`

**Usage**:
```python
GOOGLE_DRIVE_URLS = {
    "embeddings.h5": "https://drive.google.com/uc?id=YOUR_FILE_ID&export=download"
}
```

### Option 3: AWS S3 (Production)

**Best for**: Production deployments

**Steps**:
1. Create S3 bucket
2. Upload files
3. Make public or use signed URLs
4. Use boto3 to download

**Usage**:
```python
import boto3

s3 = boto3.client('s3')
s3.download_file('your-bucket', 'embeddings.h5', 'local_embeddings.h5')
```

### Option 4: File Upload in Streamlit (User-provided data)

**Best for**: Interactive use

The app already supports file uploads:
- Users can upload their own `.npy` and `.csv` files
- Files are cached locally
- No external hosting needed

## 🔧 Implementation

### 1. Update the data downloader

Edit `data_downloader.py` with your actual URLs:

```python
# Replace these with your actual URLs
GITHUB_RELEASE_URLS = {
    "full_accession_arrays.npy": "https://github.com/khoa-yelo/unnotate/releases/download/v1.0.0/full_accession_arrays.npy",
    "full_similarity_array.npy": "https://github.com/khoa-yelo/unnotate/releases/download/v1.0.0/full_similarity_array.npy",
    "parsed_uniprot_swiss_data.csv": "https://github.com/khoa-yelo/unnotate/releases/download/v1.0.0/parsed_uniprot_swiss_data.csv"
}

GOOGLE_DRIVE_URLS = {
    "embeddings.h5": "https://drive.google.com/uc?id=YOUR_EMBEDDINGS_FILE_ID&export=download"
}
```

### 2. Add data download to main app

Add this to your `streamlit_app.py`:

```python
# In the sidebar
if st.sidebar.button("📥 Download Data"):
    from data_downloader import download_sample_data
    download_sample_data()
```

## 📋 Recommended Setup

### For Streamlit Cloud Deployment:

1. **Small files** (< 100MB): Use GitHub Releases
2. **Large files** (> 100MB): Use Google Drive or AWS S3
3. **User flexibility**: Keep file upload feature

### File Structure:
```
unnotate/
├── streamlit_app.py          # Main app
├── data_downloader.py        # Data download utility
├── requirements.txt          # Dependencies
├── .streamlit/config.toml    # Streamlit config
├── .gitignore               # Excludes large files
└── data/                    # Local data cache
    ├── downloaded_files/
    └── uploaded_files/
```

## 🚀 Quick Start

1. **Upload small files to GitHub Releases**
2. **Upload large files to Google Drive**
3. **Update URLs in `data_downloader.py`**
4. **Deploy to Streamlit Cloud**
5. **Users can download data or upload their own**

## 💡 Tips

- **Cache data**: Use `@st.cache_data` for downloaded files
- **Progress indicators**: Show download progress
- **Error handling**: Graceful fallbacks if files unavailable
- **File validation**: Check file integrity after download
- **Compression**: Consider compressing large files

## 🔒 Security Notes

- **Public files**: Only upload non-sensitive data
- **Access control**: Use signed URLs for private data
- **Rate limiting**: Be mindful of download limits
- **Costs**: Monitor cloud storage costs 