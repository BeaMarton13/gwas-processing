"""
Main entry point for downloading and filtering GWAS summary statistics by disease.

Queries the GWAS Catalog spreadsheet (data.xlsx) for studies matching a disease keyword,
downloads each study's summary statistics from the EBI FTP server, filters SNPs by
p-value threshold, and saves the results as JSON files in the specified processed folder.

Usage:
    python main.py <disease> <processed_folder> <p_value_limit> <whole_word_match>

Example:
    python main.py dementia processed_0_01 0.01 True
"""
from src.utils.downloader import download_file
from src.utils.excel_handler import ExcelHandler
from src.utils.plotting import plot_keyword_bubble_chart
from src.utils.process_file import read_tsvgz_file

from pathlib import Path
import sys
import pandas as pd
import os
import json
import argparse


def download_diseases(excel_handler, disease_col, disease, disease_w_spaces=None, whole_word_match=True):
    """Download all GWAS studies matching a disease keyword from the EBI FTP server.

    Queries the GWAS Catalog spreadsheet for studies where the disease column contains
    the specified keyword, writes a manifest file listing all matching GCST IDs and
    their reported traits, then downloads each study's summary statistics file.

    Args:
        excel_handler: ExcelHandler instance wrapping the GWAS Catalog spreadsheet.
        disease_col: Column name to search (e.g., "Reported trait").
        disease: Disease keyword for filename/directory naming (underscores allowed).
        disease_w_spaces: Disease keyword for search (with spaces). If None, uses `disease`.
        whole_word_match: If True, uses regex word boundary matching; otherwise substring match.

    Outputs:
        - ./<disease>_diseases_w_id.txt — manifest of GCST IDs and trait descriptions
        - ./data/<disease>/GCST*.tsv.gz — downloaded GWAS summary statistics files
    """
    if disease_w_spaces is None:
        disease_w_spaces = disease
    study_df = excel_handler.get_disease_column(disease_col, disease_w_spaces, whole_word_match=whole_word_match)
    ids = [x for x in study_df["Study Accession"]]

    print("====================")
    print(study_df[disease_col])
    print("====================")
    print(ids)

    filename_w_path = os.path.abspath(os.path.join(sys.path[0], f"./{disease}_diseases_w_id.txt"))
    with open(filename_w_path, 'w') as f:
        # for df in study_df:
        #     print(df)
        study = study_df["Study Accession"]
        diseases = study_df[disease_col]
        for s, d in zip(study, diseases):
            f.write(f"{s} - {d}\n")

    for gcst_id in ids:
        output_path = Path(sys.path[0]) / "data" / f"{disease}" / f"{gcst_id}.tsv.gz"
        print(f"Output path: {output_path}")
        download_file(gcst_id, output_path)


def rank_keywords(excel_handler, disease_col):
    """Generate a bubble chart of keyword frequencies from reported trait text.

    Tokenizes and counts all non-stopword keywords appearing in the disease column,
    then visualizes the top N keywords as a bubble chart with size proportional to frequency.

    Args:
        excel_handler: ExcelHandler instance wrapping the GWAS Catalog spreadsheet.
        disease_col: Column name to analyze (e.g., "Reported trait").

    Returns:
        pd.DataFrame: Keyword ranking table with columns ['Keyword', 'Count'].
    """
    # Example: Get keyword ranking
    keyword_ranking = excel_handler.get_keyword_ranking(
        excel_handler.df[disease_col].dropna().astype(str)
    )
    print(keyword_ranking) 
    plot_keyword_bubble_chart(keyword_ranking, top_n=60)
    return keyword_ranking


def main():
    """Parse command-line arguments, download GWAS data, and filter by p-value threshold.

    Workflow:
        1. Parse disease name, processed folder, p-value limit, and whole-word match flag.
        2. Load the GWAS Catalog Excel file.
        3. Download all GWAS studies matching the disease keyword.
        4. For each downloaded .tsv.gz file, read in chunks and filter SNPs where
           p_value < P_VALUE_LIMIT.
        5. Save filtered SNPs as .json files in data/<disease>/<processed_folder>/.

    Command-line Args:
        disease: Disease keyword (e.g., dementia, anxiety).
        processed_folder: Output subdirectory name (e.g., processed_0_01).
        p_value_limit: Numeric p-value threshold (e.g., 0.01).
        whole_word_match: "True" or "False" — whether to use whole-word regex.
    """
    parser = argparse.ArgumentParser(
        description="Process disease-related data from Excel and TSV.GZ files."
    )
    parser.add_argument(
        "disease",
        type=str,
        help="Disease name to process (e.g., dementia, Alzheimer, Parkinson)."
    )
    parser.add_argument(
        "processed_folder",
        type=str,
        help="Processed folder where the processed files are stored (e.g. processed or processed_0_01)."
    )
    parser.add_argument(
        "p_value_limit",
        type=str,
        help="P-value limit to filter significant results."
    )
    parser.add_argument(
        "whole_word_match",
        type=str,
        help="Whether to match the disease as a whole word only (True/False)."
    )
    args = parser.parse_args()

    disease = args.disease.strip()
    disease_w_spaces = disease.replace("_", " ")
    processed_folder = args.processed_folder.strip()
    whole_word_match = args.whole_word_match.strip().lower() == "true"
    P_VALUE_LIMIT = float(args.p_value_limit.strip())

    


    if not disease:
        parser.error("❌ You must provide a disease name.\nExample: python main.py dementia")

    if not processed_folder:
        parser.error("❌ You must provide a processed folder.\nExample: python main.py dementia processed_0_01 0.01")

    if not P_VALUE_LIMIT:
        parser.error("❌ You must provide a p-value limit.\nExample: python main.py dementia processed_0_01 0.01")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    folder_path = os.path.abspath(os.path.join(script_dir, f"./data/{disease}/{processed_folder}"))
    if not os.path.exists(folder_path):
        print(f"⚠️ No processed folder found for '{disease}'. Expected at:\n{folder_path}")
        os.makedirs(folder_path)

    # === Handle Excel file ===
    file_path = sys.path[0] + "/data.xlsx"
    excel_handler = ExcelHandler(file_path)

    # === Keyword Ranking ===
    disease_col = "Reported trait" 
    pd.set_option('display.max_colwidth', None)
    pd.set_option('display.max_rows', None)
    # rank_keywords(excel_handler, disease_col)

    # === Download Disease related files ===
    download_diseases(excel_handler, disease_col, disease, disease_w_spaces=disease_w_spaces, whole_word_match=whole_word_match)

    # === Process Files ===
    # P_VALUE_LIMIT = 0.01
    # P_VALUE_LIMIT = 1 * 10 ** -5
    disease_dir = sys.path[0] + f"/data/{disease}"

    if not os.path.exists(disease_dir):
        print(f"⚠️ No data directory found for '{disease}'. Did the download step succeed?")
        sys.exit(1)

    for file in os.listdir(disease_dir):
        sample_file = os.path.join(disease_dir, file)
        if os.path.isfile(sample_file):
            print(f"Processing: {sample_file}")

            json_file = os.path.join(
                os.path.dirname(sample_file),
                folder_path,
                os.path.basename(sample_file).replace("tsv.gz", "json")
            )
            Path(json_file).parent.mkdir(parents=True, exist_ok=True)

            if not os.path.isfile(json_file) and not json_file.startswith(".") and json_file.endswith(".json"):
                selected = []
                for chunk in read_tsvgz_file(sample_file):
                    for d in chunk:
                        if float(d['p_value']) < P_VALUE_LIMIT:
                            selected.append(d)

                with open(json_file, 'w') as f:
                    json.dump(selected, f, indent=4)
                
                print(f"✅ Found {len(selected)} significant results in {sample_file}")
            else:
                print(f"ℹ️ {json_file} already exists!")

    return 


if __name__ == "__main__":
    main()