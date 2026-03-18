# import json
# import os

# disease = "dementia"

# script_dir = os.path.dirname(os.path.abspath(__file__))
# folder_path = os.path.abspath(os.path.join(script_dir, f"../../data/{disease}/processed"))

# for json_file in os.listdir(folder_path):
#     if json_file.endswith(".json"):
#         json_file_w_path = os.path.abspath(os.path.join(folder_path, json_file))
#         with open(json_file_w_path, 'r') as f:
#             json_data = json.load(f)

#         data = []
#         for d in json_data:
#             data.append({'chrom': d['chromosome'], 'start': int(d['base_pair_location']) - 1, 'end': int(d['base_pair_location'])})


#         bed_file = json_file_w_path.split('.')[0] + '.bed'
#         with open(bed_file, "w") as bed:
#             for entry in data:
#                 chrom = entry["chrom"]
#                 start = entry["start"]
#                 end = entry["end"]
#                 name = entry.get("name", ".")  # default "." if no name
#                 bed.write(f"{chrom}\t{start}\t{end}\t{name}\n")

"""
Convert filtered GWAS JSON files to BED format for genomic intersection.

Reads filtered JSON files (output from main.py) containing SNP positions and
p-values, converts them to 5-column BED format (chromosome, start, end, p_value,
name), and saves as .bed files in the same directory.

Usage:
    python convert_json_to_bed.py <disease> <processed_folder>

Example:
    python convert_json_to_bed.py dementia processed_0_01
"""
import json
import os
import argparse

def main():
    """Convert all JSON files in a processed folder to BED format.

    Parses command-line arguments for disease and processed folder, locates
    all .json files in data/<disease>/<processed_folder>/, and converts each
    to a .bed file with columns:
        chromosome | start (0-based) | end (1-based) | p_value | name

    Command-line Args:
        disease: Disease keyword (e.g., dementia, anxiety).
        processed_folder: Subdirectory name (e.g., processed_0_01).

    Outputs:
        <GCST_ID>.bed files in the same directory as the input JSON files.
    """
    parser = argparse.ArgumentParser(
        description="Convert processed JSON disease files into BED format."
    )
    parser.add_argument(
        "disease",
        type=str,
        help="Disease name to process (e.g., dementia, Alzheimer, Parkinson)."
    )
    parser.add_argument(
        "processed_folder",
        type=str,
        help="Processed folder name to use (e.g. processed or processed_0_01)."
    )
    args = parser.parse_args()
    disease = args.disease.strip()
    processed_folder = args.processed_folder.strip()

    if not disease:
        parser.error("❌ You must provide a disease name.\nExample: python convert_json_to_bed.py dementia")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    folder_path = os.path.abspath(os.path.join(script_dir, f"../../data/{disease}/{processed_folder}"))
    # folder_path = os.path.abspath(os.path.join(script_dir, f"../../data/{disease}/processed_0_01"))

    if not os.path.exists(folder_path):
        print(f"⚠️ No processed folder found for '{disease}'. Expected at:\n{folder_path}")
        return

    for json_file in os.listdir(folder_path):
        if json_file.endswith(".json") and not json_file.startswith("."):
            json_file_w_path = os.path.abspath(os.path.join(folder_path, json_file))

            with open(json_file_w_path, 'r') as f:
                json_data = json.load(f)

            data = []
            for d in json_data:
                data.append({
                    'chrom': d['chromosome'],
                    'start': int(d['base_pair_location']) - 1,
                    'end': int(d['base_pair_location']),
                    'p_value': float(d['p_value'])
                })

            bed_file = json_file_w_path.split('.')[0] + '.bed'
            with open(bed_file, "w") as bed:
                for entry in data:
                    chrom = entry["chrom"]
                    start = entry["start"]
                    end = entry["end"]
                    p_value = entry.get("p_value", ".")
                    name = entry.get("name", ".")  # default "." if no name
                    bed.write(f"{chrom}\t{start}\t{end}\t{p_value}\t{name}\n")

            print(f"✅ Created BED file: {bed_file}")


if __name__ == "__main__":
    main()