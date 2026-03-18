"""
Download GWAS summary statistics files from the NHGRI-EBI GWAS Catalog FTP server.

Constructs EBI FTP URLs from GCST accession IDs (computing the 1000-range bucket for
the directory structure) and downloads files with resumable streaming and progress bars.
"""
import requests
from pathlib import Path
import sys
from tqdm import tqdm


def build_gwas_link(gcst_id):
    """Build the EBI FTP URL for a given GCST accession ID.

    The NHGRI-EBI GWAS Catalog FTP structure groups GCST IDs into 1000-range buckets.
    For example, GCST90473238 maps to the directory GCST90473001-GCST90474000.

    Args:
        gcst_id: A string representing the GCST ID (e.g., 'GCST90473238').

    Returns:
        str: The full URL to the summary statistics .tsv.gz file.

    Raises:
        ValueError: If gcst_id does not start with 'GCST'.
    """
    if not gcst_id.startswith('GCST'):
        raise ValueError("Invalid GCST ID format. Must start with 'GCST'.")
        
    # Extract the numeric part of the ID
    numeric_id = int(gcst_id[4:])
    
    # Calculate the range for the directory name
    start_range = (numeric_id // 1000) * 1000 + 1
    end_range = start_range + 999
    
    # Format the range and the ID for the URL
    range_str = f"GCST{start_range:08}-GCST{end_range:08}"
    
    # Construct the full URL using an f-string
    base_url = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
    full_url = f"{base_url}/{range_str}/{gcst_id}/"
    filename = f"{gcst_id}.tsv.gz"

    download_url = f"{full_url}{filename}"
    return download_url

def download_file(id, local_filepath):
    """Download a GWAS summary statistics file from the EBI FTP server.

    Constructs the download URL using build_gwas_link(), ensures the parent directories
    exist, and downloads the file with streaming and a progress bar. Skips download if
    the file already exists at the target path.

    Args:
        id: GCST accession ID (e.g., 'GCST90473238').
        local_filepath: Full local path where the file should be saved (e.g., Path object or str).

    Outputs:
        Writes the downloaded .tsv.gz file to local_filepath if not already present.
        Prints status messages to stdout.
    """
    download_url = build_gwas_link(id)
    print(f"Downloading from: {download_url}")
    # Ensure parent directory exists
    Path(local_filepath).parent.mkdir(parents=True, exist_ok=True)
    Path(local_filepath).parent.parent.mkdir(parents=True, exist_ok=True)

    if not Path(local_filepath).is_file():
        try:
            with requests.get(download_url, stream=True) as response:
                response.raise_for_status()
                total_size = int(response.headers.get('content-length', 0))
                chunk_size = 1024
                with tqdm(total=total_size, unit='B', unit_scale=True, desc="Downloading") as pbar:
                    with open(local_filepath, 'wb') as f:
                        for chunk in response.iter_content(chunk_size=chunk_size):
                            if chunk:
                                f.write(chunk)
                                pbar.update(len(chunk))
            print(f"File successfully downloaded and saved as: {local_filepath}")
        except requests.exceptions.HTTPError as errh:
            print(f"HTTP Error: {errh}")
        except requests.exceptions.ConnectionError as errc:
            print(f"Error Connecting: {errc}")
        except requests.exceptions.Timeout as errt:
            print(f"Timeout Error: {errt}")
        except requests.exceptions.RequestException as err:
            print(f"An error occurred: {err}")
    else:
        print(f"File already exists at: {local_filepath}")


if __name__ == "__main__":
    gcst_id = "GCST90473238"
    gcst_id = "GCST90473304"
    
    output_path = Path(sys.path[0]).parent.parent / "data" / f"{gcst_id}.tsv.gz"
    print(f"Output path: {output_path}")
    download_file(gcst_id, output_path)
