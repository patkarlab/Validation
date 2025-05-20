import json
import pandas as pd
import sys
import os

args = sys.argv

input_json_file = args[1]
output_csv_file_details = args[2]

try:
    # Check if the input JSON file is empty
    if os.stat(input_json_file).st_size == 0:
        print(f"Input JSON file {input_json_file} is empty. Creating an empty CSV file: {output_csv_file_details}")
        pd.DataFrame().to_csv(output_csv_file_details, index=False)
        sys.exit(0)

    # Read the JSON file and replace problematic values
    with open(input_json_file, "r", encoding="utf-8") as file:
        raw_content = file.read()
        # Replace invalid JSON values
        raw_content = raw_content.replace("-nan", "null")

    # Load the cleaned JSON data
    data = json.loads(raw_content)
except UnicodeDecodeError:
    with open(input_json_file, "r", encoding="latin1") as file:
        raw_content = file.read()
        raw_content = raw_content.replace("-nan", "null")
    data = json.loads(raw_content)

# Check if the "details" key exists and is non-empty
if "details" not in data or not data["details"]:
    print(f"'Details' key is empty in JSON file {input_json_file}. Creating an empty CSV file: {output_csv_file_details}")
    pd.DataFrame().to_csv(output_csv_file_details, index=False)
    sys.exit(0)

# Proceed if "details" contains data
details_df = pd.DataFrame(data["details"])

if "VAF" in details_df.columns:
    details_df["VAF"] *= 100

details_df.to_csv(output_csv_file_details, index=False)
print(f"Details saved to: {output_csv_file_details}")
