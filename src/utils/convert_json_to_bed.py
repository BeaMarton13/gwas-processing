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

import json
import os
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Convert processed JSON disease files into BED format."
    )
    parser.add_argument(
        "disease",
        type=str,
        help="Disease name to process (e.g., dementia, Alzheimer, Parkinson)."
    )
    args = parser.parse_args()
    disease = args.disease.strip()

    if not disease:
        parser.error("❌ You must provide a disease name.\nExample: python convert_json_to_bed.py dementia")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    folder_path = os.path.abspath(os.path.join(script_dir, f"../../data/{disease}/processed"))

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
                    'end': int(d['base_pair_location'])
                })

            bed_file = json_file_w_path.split('.')[0] + '.bed'
            with open(bed_file, "w") as bed:
                for entry in data:
                    chrom = entry["chrom"]
                    start = entry["start"]
                    end = entry["end"]
                    name = entry.get("name", ".")  # default "." if no name
                    bed.write(f"{chrom}\t{start}\t{end}\t{name}\n")

            print(f"✅ Created BED file: {bed_file}")


if __name__ == "__main__":
    main()