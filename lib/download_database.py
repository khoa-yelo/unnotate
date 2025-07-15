import gdown
import os

def download_gdrive_folder(dest_path: str = "./database"):
    """
    Download an entire folder from Google Drive using gdown.
    Args:
        dest_path (str): Local directory to save the downloaded folder.
    """
    url = f"https://drive.google.com/drive/folders/1_kLGss-g__JIwdVzNriB18T9pUsklmxS"
    os.makedirs(dest_path, exist_ok=True)
    try:
        gdown.download_folder(url, output=dest_path, quiet=False, use_cookies=False)
        print(f"Download complete. Folder saved to: {dest_path}")
    except Exception as e:
        print(f"Error downloading folder: {e}")

if __name__ == "__main__":
    download_gdrive_folder()
