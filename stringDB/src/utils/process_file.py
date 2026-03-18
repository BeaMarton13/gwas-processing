"""
Read and parse GWAS summary statistics files and GTF genome annotations.

Provides chunked readers for large gzipped TSV files (memory-efficient for 3+ GB
GWAS summary stats) and utilities for sorting SNPs by p-value.
"""
import sys
import json
import gzip
import pandas as pd

from gtfparse import read_gtf
from pathlib import Path

P_VALUE_LIMIT = 0.01
# P_VALUE_LIMIT = 1 * 10 ** -5


def read_tsvgz_file(file_path, chunk_size=1000):
    """Read a gzipped TSV file in chunks and yield rows as lists of dictionaries.

    Memory-efficient generator for large GWAS summary statistics files. Reads
    the header line to extract column names, then yields rows in batches of
    chunk_size dictionaries.

    Args:
        file_path: Path to a .tsv.gz file.
        chunk_size: Number of rows to yield per chunk (default: 1000).

    Yields:
        list[dict]: Chunks of rows, where each row is a dict mapping column names to values.
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
    """Read a GTF file and return its contents as JSON.

    Uses the gtfparse library to parse the GTF format and converts the
    resulting DataFrame to JSON.

    Args:
        file_path: Path to a GTF file (can be gzipped).

    Returns:
        str: JSON-serialized GTF data, or None if an error occurs.
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
    """Sort a list of SNP dictionaries by p-value in ascending order.

    Args:
        data: List of dictionaries, each containing a 'p_value' key.

    Returns:
        list[dict]: The sorted list (smallest p-values first).
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
