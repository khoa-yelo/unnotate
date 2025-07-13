# Streamlit Cloud Deployment Guide

This guide will help you deploy the Protein Accession Visualizer to Streamlit Cloud with automatic data loading from GitHub releases.

## Prerequisites

1. ✅ GitHub repository with code pushed
2. ✅ GitHub release with data files (`protein_data_v1.0.zip`)
3. ✅ Streamlit Cloud account

## Step-by-Step Deployment

### 1. Go to Streamlit Cloud
Visit [share.streamlit.io](https://share.streamlit.io) and sign in with your GitHub account.

### 2. Create New App
1. Click "New app"
2. Select your repository: `khoa-yelo/unnotate`
3. Set the main file path to: `streamlit_app.py`
4. Click "Deploy"

### 3. Configure Secrets (Optional)
The app will automatically download data from your GitHub release. If you want to customize the repository settings:

1. Go to your deployed app
2. Click the hamburger menu (☰) → "Settings"
3. Scroll down to "Secrets"
4. Add these secrets:
   ```toml
   GITHUB_REPO_OWNER = "khoa-yelo"
   GITHUB_REPO_NAME = "unnotate"
   GITHUB_RELEASE_TAG = "latest"
   ```

### 4. How It Works

When the app starts on Streamlit Cloud:

1. **Checks for data files** in the repository
2. **If missing**: Automatically downloads from your GitHub release
3. **Extracts the zip file** and loads the data
4. **Displays the visualization** with all protein data

## Troubleshooting

### Deployment Fails
- **Check repository**: Ensure all code is pushed to GitHub
- **Check main file**: Verify `streamlit_app.py` exists in the root
- **Check requirements**: Ensure `requirements.txt` is present

### Data Download Fails
- **Check GitHub release**: Ensure `protein_data_v1.0.zip` exists
- **Check repository name**: Verify it matches your GitHub repo
- **Check internet**: Streamlit Cloud needs internet access

### App Loads but No Data
- **Check logs**: Look for download errors in the Streamlit logs
- **Manual download**: Use the download interface in the app
- **Verify release**: Check that the zip file contains all required files

## Expected Behavior

1. **First deployment**: App will download ~135MB of data (takes 1-2 minutes)
2. **Subsequent visits**: Data is cached, loads quickly
3. **Data source**: Always uses the latest release from your GitHub repository

## File Structure for Cloud

```
unnotate/
├── streamlit_app.py          # Main app (auto-downloads data)
├── requirements.txt          # Dependencies
├── .streamlit/
│   ├── config.toml          # Streamlit config
│   └── secrets.toml         # GitHub repo settings
├── README.md                # Documentation
└── (data files downloaded at runtime)
```

## Benefits of This Approach

- ✅ **No large files in repository** (keeps it lightweight)
- ✅ **Automatic updates** (new releases automatically available)
- ✅ **Scalable** (works for any size dataset)
- ✅ **User-friendly** (no manual data management needed)

## Support

If you encounter issues:
1. Check the Streamlit Cloud logs
2. Verify your GitHub release exists
3. Test the download manually using the app interface 