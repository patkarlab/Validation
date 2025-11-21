#!/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Take the fltr3.json.csv and get_itd outfiles as input files, find common variant between flt3r and get_itd, add to combined output excel file along with the flt3r annovar_annotated file, varscan output and coverage file")
parser.add_argument("-s", "--sample_id", required=True, help="sample_id to be used as prefix")
parser.add_argument("-f", "--flt3r_json", required=True, help="Path to flt3r output")
parser.add_argument("-m", "--merged_out", required=True, help="Path to merged varscan & flt3r annovar-annotated output")
parser.add_argument("-g", "--getitd_out", required=True, help="Path to get_itd output")
parser.add_argument("-c", "--cov_out", required=True, help="Path to coverage file")
parser.add_argument("-o", "--output", required=True, help="Path to out excel file")

args = parser.parse_args()

#read in merged varscan and filt3r out, and coverage out file
merged_sheet1 = pd.read_excel(args.merged_out, sheet_name="varscan")
merged_sheet2 = pd.read_excel(args.merged_out, sheet_name="filt3r")
cov_bed = pd.read_csv(args.cov_out, sep = "\t", names=['Chrom', 'Start', 'Stop', 'Region', 'Coverage'])

#read in input files
df1 = pd.read_csv(args.flt3r_json)
df_flt3r_in = df1.copy()
df2 = pd.read_csv(args.getitd_out, sep="\t")
df_getITD_in = df2.copy()

#filter df1
#df1 = df1.loc[df1['is_wt_duplication'] == True]
df1 = df1.loc[df1['size'] >= 5 ]
df1 = df1.loc[df1['occurrence'] >= 2 ]

#filter df2
df2 = df2.loc[df2['length'] >= 5 ]
df2 = df2.loc[df2['counts'] >= 2 ]

#add tool name suffix to each column
df1 = df1.add_suffix('_FiLT3r')
df2 = df2.add_suffix('_getITD')

#merge on common values in size and length of each
merged_df = pd.merge(df1, df2, left_on='size_FiLT3r', right_on='length_getITD', how='inner')

#add a new column for ALT counts from getITD
merged_df['Ref_Counts_getITD'] = merged_df['coverage_getITD'] - merged_df['counts_getITD']

merged_df = merged_df[['sample_getITD','length_getITD','counts_getITD', 'Ref_Counts_getITD', 'seq_getITD','start_chr13_bp_getITD','end_chr13_bp_getITD', 'vaf_getITD','size_FiLT3r','occurrence_FiLT3r' , 'wt_coverage_FiLT3r', 'VAF_FiLT3r','sequence_FiLT3r','reference_pos_FiLT3r','is_wt_duplication_FiLT3r' ]]

#group rows which have same length & sequence in getITD and keep rows with top2 vaf from flt3r
merged_df = merged_df.groupby(['length_getITD','seq_getITD']).apply(lambda x: x.sort_values('occurrence_FiLT3r', ascending=False).head(2))

merged_df.rename(columns={'length_getITD':'ITD_length_getITD', 'counts_getITD':'Alt_counts_getITD', 'start_chr13_bp_getITD':'start_getITD', 'end_chr13_bp_getITD':'end_getITD','vaf_getITD':'VAF_getITD(%)', 'size_FiLT3r':'ITD_length_FiLT3r', 'occurrence_FiLT3r':'Alt_Counts_FiLT3r', 'wt_coverage_FiLT3r':'Ref_Counts_FiLT3r', 'VAF_FiLT3r':'VAF_FiLT3r(%)', 'reference_pos_FiLT3r':'chr_start_strand_FiLT3r'   }, inplace=True)

#print outfile
with pd.ExcelWriter(args.output) as writer:
	merged_df.to_excel(writer, sheet_name = args.sample_id + "_common", index=False)
	df_getITD_in.to_excel(writer, sheet_name= args.sample_id + "_getITD", index=False)
	df_flt3r_in.to_excel(writer, sheet_name= args.sample_id + "_FiLT3r_json", index=False)
	merged_sheet2.to_excel(writer, sheet_name= args.sample_id +"_Filt3r_anno", index=False)
	merged_sheet1.to_excel(writer, sheet_name= args.sample_id + "_varscan", index=False)
	cov_bed.to_excel(writer, sheet_name="coverage", index=False)



