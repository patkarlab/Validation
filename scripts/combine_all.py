#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Take the combined output, fltr3 and get_itd outfiles as input files, find common variant between flt3r and get_itd, add to combined output excel file")

#parser.add_argument("--flt3r_out", required=True, help="Path to flt3r output")
parser.add_argument("--getitd_out", required=True, help="Path to get_itd output")
parser.add_argument("--merged_out", required=True, help="Path to merged varscan and flt3r output")
parser.add_argument("--cov_out", required=True, help="Path to coverage file")
parser.add_argument("--output", required=True, help="Path to out excel file")

args = parser.parse_args()

#read in input files
merged_sheet1 = pd.read_excel(args.merged_out, sheet_name="varscan")
df1 = pd.read_excel(args.merged_out, sheet_name="filt3r")
df2 = pd.read_csv(args.getitd_out, sep = "\t" )
cov_bed = pd.read_csv(args.cov_out, sep = "\t", names=['Chrom', 'Start', 'Stop', 'Region', 'Coverage'])
#df1 = pd.read_csv(args.flt3r_out)

#add tool name suffix to each column 
df1 = df1.add_suffix('_flt3r')
df2 = df2.add_suffix('_getitd')

#join the 2 dataframes on the basis of key, drop the key column
df3 = df1.assign(key=0).merge(df2.assign(key=0), how='left', on = 'key').drop(columns=['key'])

#function to find the longest common substring(min length=5) between 2 sequences
def common_substring(seq1, seq2):
    min_len = 5
    max_len = 0
    match_seq = ""

    if len(seq1) == len(seq2):
        for i in range (len (seq1)):
            for j in range (1 + min_len, len(seq1) +1):
                substr = seq1[i:j]
                if substr in seq2 and len(substr) >= min_len:
                    if len(substr ) > max_len:
                        max_len = len(substr)
                        match_seq = substr


    if len(match_seq) > min_len:
        return match_seq
    else:
        return False

try:
    required_cols = ['sequence_flt3r', 'seq_getitd']

    if not df3.empty:
        # Check if required columns exist and are not fully empty (i.e., not all NaN)
        if all(col in df3.columns for col in required_cols) and \
           not df3[required_cols].isna().all().any():

            #apply the common_substring function on the columns containing itd sequence
            df3['matching_itd'] = df3.apply(lambda row: common_substring(row["sequence_flt3r"], row["seq_getitd"]), axis = 1)

        else:
            print("Skipping processing: Required columns exist but contain only NaNs.")
            df3['matching_itd'] = False # optional placeholder
    else:
        print("Skipping processing: Merged DataFrame is completely empty.")
        df3['matching_itd'] = False   # optional placeholder

except Exception as e:
    print(f"Error during processing: {e}")
    if 'matching_itd' not in df3.columns:
        df3['matching_itd'] = False   # avoid breaking Excel export


#remove the column containing the common seq between flt3r and get_itd itds and insert the column at index 0
common_itd = df3.pop('matching_itd')
df3.insert(0, 'matching_itd', common_itd)

#extract rows where common variant identified
match_df = df3[df3['matching_itd'] != False ]

#print merged output
with pd.ExcelWriter(args.output) as writer:
    merged_sheet1.to_excel(writer, sheet_name="varscan", index=False)
    df1.to_excel(writer, sheet_name="flt3r", index=False)
    df2.to_excel(writer, sheet_name="get_itd", index=False)
    match_df.to_excel(writer, sheet_name="common_itds", index=False)
    cov_bed.to_excel(writer, sheet_name="coverage", index=False)

