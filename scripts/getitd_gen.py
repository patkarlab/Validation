import pandas as pd
import argparse

columns = [
    "sample",
    "length",
    "start",
    "vaf",
    "ar",
    "coverage",
    "counts",
    "trailing",
    "seq",
    "sense",
    "external_bp",
    "domains",
    "start_chr13_bp",
    "start_transcript_bp",
    "start_protein_as",
    "end_chr13_bp",
    "end_transcript_bp",
    "end_protein_as",
    "insertion_site_chr13_bp",
    "insertion_site_transcript_bp",
    "insertion_site_protein_as",
    "insertion_site_domain",
    "file"
]

parser = argparse.ArgumentParser(description='Create an empty table with specified columns.')
parser.add_argument('-o', '--output', default='empty_table.tsv', help='Output file name (default: empty_table.tsv)')
args = parser.parse_args()

df = pd.DataFrame(columns=columns)
df.to_csv(args.output, sep="\t", index=False)
print(f"Empty table saved to {args.output}")
