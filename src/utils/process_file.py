import sys
import json
import gzip
import pandas as pd

from gtfparse import read_gtf
from pathlib import Path

P_VALUE_LIMIT = 1 * 10 ** -5


# def read_tsvgz_file(file_path):
#     """
#     Reads a gzipped TSV file and returns its content as a list of dictionaries.
    
#     Args:
#         file_path (str): The path to the gzipped TSV file.
#     """
#     data = []
#     with gzip.open(file_path, 'rb') as file:
#         headers = [h.decode() for h in file.readline().strip().split(b'\t')]
#         idx = 0
#         for line in file:
#             values = [v.decode() for v in line.strip().split(b'\t')]
#             entry = dict(zip(headers, values))
#             data.append(entry)
#             # idx += 1
#             # if idx == 50:
#             #     break
#     return data


def read_tsvgz_file(file_path, chunk_size=1000):
    """
    Reads a gzipped TSV file and yields chunks of data as lists of dictionaries.
    
    Args:
        file_path (str): The path to the gzipped TSV file.
        chunk_size (int): Number of lines to process at once.
    """
    with gzip.open(file_path, 'rb') as file:
        headers = [h.decode() for h in file.readline().strip().split(b'\t')]
        chunk = []
        
        for line in file:
            values = [v.decode() for v in line.strip().split(b'\t')]
            entry = dict(zip(headers, values))
            chunk.append(entry)
            
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        
        if chunk:  # Don't forget the last partial chunk
            yield chunk



def read_gtf_file(file_path):
    """
    Reads a GTF file and returns a pandas DataFrame.
    
    Args:
        file_path (str): The path to the GTF file.
    
    Returns:
        pandas.DataFrame: A DataFrame containing the parsed GTF data.
    """
    try:
        # read_gtf handles the parsing of the fixed columns and the attributes
        df = read_gtf(file_path)
        print(type(df))
        return df.write_json()
    except Exception as e:
        print(f"Error reading GTF file: {e}")
        return None


def sort_by_p_value(data):
    """
    Sorts a list of dictionaries by the 'P-value' key in ascending order.
    
    Args:
        data (list): A list of dictionaries, each containing a 'P-value' key.
        
    Returns:
        list: The sorted list of dictionaries.
    """
    return sorted(data, key=lambda x: float(x['p_value']))


if __name__ == "__main__":
    sample_file = sys.path[0] + "/../../data/mental/GCST90473237.tsv.gz"
    data = read_tsvgz_file(sample_file)

    json_file = "/".join(sample_file.split("/")[:-1]) + "/processed/" + sample_file.split("/")[-1].replace("tsv.gz", "json")
    Path(json_file).parent.mkdir(parents=True, exist_ok=True)

    selected = []
    for d in data:
        print(d['p_value'])
        if float(d['p_value']) < P_VALUE_LIMIT:
            selected.append(d)

    with open(json_file, 'w') as f:
        json.dump(selected, f, indent=4) 
