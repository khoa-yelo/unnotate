"""
Data downloader utility for the protein accession visualizer
"""
import streamlit as st
import requests
import zipfile
import os
import tempfile
import shutil
from pathlib import Path

def download_from_github_release(repo_owner, repo_name, release_tag="latest"):
    """
    Download data files from a GitHub release
    
    Args:
        repo_owner: GitHub username/organization
        repo_name: Repository name
        release_tag: Release tag (default: "latest")
    """
    
    if release_tag == "latest":
        # Get latest release info
        api_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/releases/latest"
    else:
        # Get specific release info
        api_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/releases/tags/{release_tag}"
    
    try:
        response = requests.get(api_url)
        response.raise_for_status()
        release_data = response.json()
        
        # Find the zip file asset
        zip_asset = None
        for asset in release_data.get('assets', []):
            if asset['name'].endswith('.zip'):
                zip_asset = asset
                break
        
        if not zip_asset:
            st.error("No zip file found in the release")
            return False
        
        # Download the zip file
        zip_url = zip_asset['browser_download_url']
        st.info(f"Downloading {zip_asset['name']} from GitHub release...")
        
        with st.spinner("Downloading data files..."):
            zip_response = requests.get(zip_url, stream=True)
            zip_response.raise_for_status()
            
            # Save to temporary file
            with tempfile.NamedTemporaryFile(delete=False, suffix='.zip') as tmp_file:
                for chunk in zip_response.iter_content(chunk_size=8192):
                    tmp_file.write(chunk)
                tmp_path = tmp_file.name
            
            # Extract to current directory
            with zipfile.ZipFile(tmp_path, 'r') as zip_ref:
                zip_ref.extractall('.')
            
            # Clean up temporary file
            os.unlink(tmp_path)
        
        st.success("✅ Data files downloaded and extracted successfully!")
        return True
        
    except requests.exceptions.RequestException as e:
        st.error(f"Failed to download from GitHub: {e}")
        return False
    except Exception as e:
        st.error(f"Error processing download: {e}")
        return False

def upload_data_files():
    """Upload data files using Streamlit file uploader"""
    st.subheader("Upload Data Files")
    
    uploaded_files = st.file_uploader(
        "Upload your data files",
        type=['npy', 'csv'],
        accept_multiple_files=True,
        help="Upload .npy files for accession arrays and similarity data, and .csv for protein metadata"
    )
    
    if uploaded_files:
        for uploaded_file in uploaded_files:
            # Save the uploaded file
            with open(uploaded_file.name, "wb") as f:
                f.write(uploaded_file.getbuffer())
            st.success(f"✅ {uploaded_file.name} uploaded successfully!")
        
        return True
    return False

def main():
    st.title("Data File Manager")
    
    st.write("""
    This tool helps you download or upload the data files needed for the protein accession visualizer.
    
    **Required files:**
    - `full_accession_arrays.npy` - Protein accession arrays
    - `full_similarity_array.npy` - Similarity scores  
    - `parsed_uniprot_swiss_data.csv` - Protein metadata
    """)
    
    # Check if files already exist
    existing_files = []
    required_files = [
        "full_accession_arrays.npy",
        "full_similarity_array.npy", 
        "parsed_uniprot_swiss_data.csv"
    ]
    
    for file in required_files:
        if os.path.exists(file):
            existing_files.append(file)
    
    if existing_files:
        st.success(f"✅ Found existing files: {', '.join(existing_files)}")
    
    # Download from GitHub release
    st.subheader("Download from GitHub Release")
    
    col1, col2 = st.columns(2)
    
    with col1:
        repo_owner = st.text_input("Repository Owner", value="your-username")
        repo_name = st.text_input("Repository Name", value="your-repo-name")
    
    with col2:
        release_tag = st.text_input("Release Tag", value="latest", 
                                   help="Use 'latest' for the most recent release, or specify a tag like 'v1.0.0'")
        
        if st.button("Download from GitHub"):
            if repo_owner and repo_name:
                success = download_from_github_release(repo_owner, repo_name, release_tag)
                if success:
                    st.rerun()  # Refresh to show updated file status
            else:
                st.error("Please enter repository owner and name")
    
    # Upload files
    st.subheader("Upload Files Manually")
    if st.button("Upload Files"):
        upload_data_files()

if __name__ == "__main__":
    main() 