#!/usr/bin/env python3
"""
Simple script to create GitHub releases for data files
"""
import os
import sys
import subprocess
import zipfile
from pathlib import Path

def create_data_zip():
    """Create a zip file with all necessary data files including CSV"""
    data_files = [
        "full_accession_arrays.npy",
        "full_similarity_array.npy", 
        "viral_accession_arrays.npy",
        "viral_similarity_array.npy",
        "bacterial_accession_arrays.npy", 
        "bacterial_similarity_array.npy"
    ]
    
    # Add CSV file if it exists
    csv_path = "../data/unnotate/parsed_uniprot_swiss_data.csv"
    if os.path.exists(csv_path):
        data_files.append(csv_path)
    
    zip_name = "protein_data_v1.0.zip"
    
    with zipfile.ZipFile(zip_name, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for file_path in data_files:
            if os.path.exists(file_path):
                print(f"Adding {file_path} to zip...")
                zipf.write(file_path, os.path.basename(file_path))
            else:
                print(f"Warning: {file_path} not found, skipping...")
    
    print(f"Created {zip_name}")
    return zip_name



def create_github_release():
    """Create a GitHub release with the data zip file"""
    # Create the zip file with all data
    zip_name = create_data_zip()
    
    # Get current version from git tags or use default
    try:
        result = subprocess.run(['git', 'describe', '--tags', '--abbrev=0'], 
                              capture_output=True, text=True, check=True)
        current_version = result.stdout.strip()
        # Increment version number
        if current_version.startswith('v'):
            version_num = current_version[1:]
            try:
                major, minor, patch = version_num.split('.')
                new_patch = int(patch) + 1
                new_version = f"v{major}.{minor}.{new_patch}"
            except:
                new_version = "v1.0.1"
        else:
            new_version = "v1.0.1"
    except subprocess.CalledProcessError:
        new_version = "v1.0.0"
    
    print(f"Creating release {new_version}...")
    
    # Create git tag
    subprocess.run(['git', 'tag', new_version], check=True)
    
    # Push tag to GitHub
    subprocess.run(['git', 'push', 'origin', new_version], check=True)
    
    # Create release using GitHub CLI (if available)
    try:
        files_to_upload = [zip_name]
        
        subprocess.run([
            'gh', 'release', 'create', new_version,
            *files_to_upload,
            '--title', f'Protein Data Release {new_version}',
            '--notes', f'''
# Protein Accession Data Release {new_version}

This release contains the protein accession data files needed for the Streamlit visualization app.

## Files included:

        ### Complete Data Package (`protein_data_v1.0.zip`)
        - `full_accession_arrays.npy` - Complete protein accession arrays
        - `full_similarity_array.npy` - Complete similarity scores
        - `viral_accession_arrays.npy` - Viral protein accession arrays  
        - `viral_similarity_array.npy` - Viral similarity scores
        - `bacterial_accession_arrays.npy` - Bacterial protein accession arrays
        - `bacterial_similarity_array.npy` - Bacterial similarity scores
        - `parsed_uniprot_swiss_data.csv` - Protein metadata from UniProt (required for full functionality)

## Usage:

1. **Download the zip file** and extract it to your project directory
2. **Install dependencies**: `pip install -r requirements.txt`
3. **Run the app**: `streamlit run streamlit_app.py`

The app will automatically detect and load these files.

## Note:
The CSV file is large (~680MB) but is required for the app to display protein information, domain classification, and other metadata. The app will not work properly without this file.
            '''.strip()
        ], check=True)
        print(f"✅ GitHub release {new_version} created successfully!")
        print(f"📦 Data files uploaded: {', '.join(files_to_upload)}")
        
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("⚠️  GitHub CLI not found or failed. Please create the release manually:")
        print(f"1. Go to your GitHub repository")
        print(f"2. Click 'Releases' → 'Create a new release'")
        print(f"3. Tag: {new_version}")
        print(f"4. Title: Protein Data Release {new_version}")
        print(f"5. Upload the files: {', '.join(files_to_upload)}")
        print(f"6. Publish the release")

if __name__ == "__main__":
    create_github_release() 