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


def download_diseases(excel_handler, disease_col, disease):
    study_df = excel_handler.get_disease_column(disease_col, disease)
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
    # Example: Get keyword ranking
    keyword_ranking = excel_handler.get_keyword_ranking(
        excel_handler.df[disease_col].dropna().astype(str)
    )
    print(keyword_ranking) 
    plot_keyword_bubble_chart(keyword_ranking, top_n=60)
    return keyword_ranking


def main():
    parser = argparse.ArgumentParser(
        description="Process disease-related data from Excel and TSV.GZ files."
    )
    parser.add_argument(
        "disease",
        type=str,
        help="Disease name to process (e.g., dementia, Alzheimer, Parkinson)."
    )
    args = parser.parse_args()

    disease = args.disease.strip()

    if not disease:
        parser.error("❌ You must provide a disease name.\nExample: python main.py dementia")

    # === Handle Excel file ===
    file_path = sys.path[0] + "/data.xlsx"
    excel_handler = ExcelHandler(file_path)

    # === Keyword Ranking ===
    disease_col = "Reported trait" 
    pd.set_option('display.max_colwidth', None)
    pd.set_option('display.max_rows', None)
    # rank_keywords(excel_handler, disease_col)

    # === Download Disease related files ===
    download_diseases(excel_handler, disease_col, disease)

    # === Process Files ===
    P_VALUE_LIMIT = 1 * 10 ** -5
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
                "processed",
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