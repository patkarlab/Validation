import json
import pandas as pd
import sys

args = sys.argv

input_json_file = args[1]
output_excel_file = args[2]

try:
    with open(input_json_file, "r", encoding="utf-8") as file:
        data = json.load(file)
except UnicodeDecodeError:
    with open(input_json_file, "r", encoding="latin1") as file:
        data = json.load(file)

details_df = pd.DataFrame(data["details"])

if "VAF" in details_df.columns:
    details_df["VAF"] *= 100

metadata = {k: v for k, v in data.items() if k != "details"}
metadata_df = pd.DataFrame(list(metadata.items()), columns=["Key", "Value"])

with pd.ExcelWriter(output_excel_file) as writer:
    details_df.to_excel(writer, sheet_name="Details", index=False)
    metadata_df.to_excel(writer, sheet_name="Metadata", index=False)

print(f"Excel file created: {output_excel_file}")
