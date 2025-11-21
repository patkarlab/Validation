#!/usr/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Take fltr3.json.csv, flt3_itd_ext and get_itd outfiles as input files, find common variant between 3 ITD callers, add to combined output excel file along with the flt3r annovar_annotated file, varscan output and coverage file")
parser.add_argument("-s", "--sample_id", required=True, help="sample_id to be used as prefix")
parser.add_argument("-j", "--flt3r_json", required=True, help="Path to flt3r output")
parser.add_argument("-f", "--flt3r_anno", required=True, help="Path to flt3r annovar-annotated output")
parser.add_argument("-v", "--varscan_out", required=True, help="Path to varscan output")
parser.add_argument("-g", "--getitd_out", required=True, help="Path to get_itd output")
parser.add_argument("-x", "--flt3_itd_ext_out", required=True, help="Path to flt3_itd_ext output")
parser.add_argument("-c", "--cov_out", required=True, help="Path to coverage file")
parser.add_argument("-o", "--output", required=True, help="Path to out excel file")

args = parser.parse_args()

#read varscan, flt3r annotated vcf and coverage files(to be merged in final excel )
try:
	varscan_df = pd.read_csv(args.varscan_out)
	if varscan_df.shape[0] == 0:
		varscan_df = pd.DataFrame(columns=['Chr','Start','End','Ref','Alt','FILTER','SOMATIC_FLAG','REF_COUNT','ALT_COUNT',
		'VAF%','Func.refGene','Gene.refGene','ExonicFunc.refGene','AAChange.refGene','Gene_full_name.refGene',
		'Function_description.refGene','Disease_description.refGene','cosmic84','PopFreqMax'])
except pd.errors.EmptyDataError:
	print(f"The file '{args.varscan_out}' is completely empty. Creating empty DataFrame with headers")
	varscan_df = pd.DataFrame(columns=['Chr','Start','End','Ref','Alt','FILTER','SOMATIC_FLAG','REF_COUNT','ALT_COUNT',
	'VAF%','Func.refGene','Gene.refGene','ExonicFunc.refGene','AAChange.refGene','Gene_full_name.refGene',
	'Function_description.refGene','Disease_description.refGene','cosmic84','PopFreqMax'])

flt3r_df = pd.read_csv(args.flt3r_anno)
cov_bed = pd.read_csv(args.cov_out, sep = "\t", names=['Chrom', 'Start', 'Stop', 'Region', 'Coverage'])

#read in flt3r json, getITD and flt3_itd_ext files (for getting common ITDs)
df1 = pd.read_csv(args.flt3r_json)
#copy the flt3r json (to be merged in final excel)
df_flt3r_in = df1.copy()

df2 = pd.read_csv(args.getitd_out, sep="\t", skiprows=1, index_col=False, names=['sample','length','start','vaf','ar','coverage','counts','trailing','seq','sense','domains','start_chr13_bp','start_transcript_bp','start_protein_as','end_chr13_bp',
'end_transcript_bp','end_protein_as','insertion_site_chr13_bp','insertion_site_transcript_bp','insertion_site_protein_as',	
'insertion_site_domain','file'])

df2 = df2.sort_values(by='vaf', ascending=False)
#copy the get_ITD file (to be merged in final excel)
df_getITD_in = df2.copy()

df3 = pd.read_csv(args.flt3_itd_ext_out, sep ="\t")
#copy the flt3_itd_ext file (to be merged in final excel)
df_flt3_itd_ext_in = df3.copy()

#finding common variants
merged_df = pd.DataFrame()
proceed_merge = True

if df1.empty or df2.empty or df3.empty:
		print("⚠️ One or more input DataFrames (df1, df2, df3) are empty. Skipping merging.")
		proceed_merge = False

else:
		# --- Step 3: Perform sequential merges ---
		try:		
			#filter df1
			df1 = df1.loc[df1['size'] >= 5 ]
			df1 = df1.loc[df1['occurrence'] >= 2 ]

			#filter df2
			df2 = df2.loc[df2['length'] >= 5 ]
			df2 = df2.loc[df2['counts'] >= 2 ]

			#filter df3
			df3 = df3.loc[df3['length'] >= 5 ]
			df3 = df3.loc[df3['Raw Variant Depth(RVD)'] >= 2 ]

			#add tool name suffix to each column
			df1 = df1.add_suffix('_FiLT3r')
			df2 = df2.add_suffix('_getITD')
			df3 = df3.add_suffix('_flt3_ext')

			#merge on common values in size and length of each
			merged_df = pd.merge(df1, df2, left_on='size_FiLT3r', right_on='length_getITD', how='inner')
			
			#add a new column for ALT counts from getITD
			merged_df['Ref_Counts_getITD'] = merged_df['coverage_getITD'] - merged_df['counts_getITD']

			merged_df = merged_df[['sample_getITD','length_getITD','counts_getITD', 'Ref_Counts_getITD', 'seq_getITD',
			'start_chr13_bp_getITD','end_chr13_bp_getITD', 'vaf_getITD','size_FiLT3r','occurrence_FiLT3r' ,
			'wt_coverage_FiLT3r', 'VAF_FiLT3r','sequence_FiLT3r','reference_pos_FiLT3r','is_wt_duplication_FiLT3r' ]]

			#group rows which have same length & sequence in getITD and keep rows with top2 vaf from flt3r
			merged_df = merged_df.groupby(['length_getITD','seq_getITD']).apply(lambda x: x.sort_values('occurrence_FiLT3r',
			ascending=False).head(2))

			#rename columns
			merged_df.rename(columns={'length_getITD':'ITD_length_getITD', 'counts_getITD':'Alt_counts_getITD',
			'start_chr13_bp_getITD':'start_getITD', 'end_chr13_bp_getITD':'end_getITD','vaf_getITD':'VAF_getITD(%)',
			'size_FiLT3r':'ITD_length_FiLT3r', 'occurrence_FiLT3r':'Alt_Counts_FiLT3r',
			'wt_coverage_FiLT3r':'Ref_Counts_FiLT3r', 'VAF_FiLT3r':'VAF_FiLT3r(%)',
			'reference_pos_FiLT3r':'chr_start_strand_FiLT3r'}, inplace=True)

			merged_df = pd.merge(merged_df, df3, left_on='ITD_length_getITD', right_on='length_flt3_ext', how='inner')

			merged_df = merged_df.groupby(['ITD_length_getITD','seq_getITD']).apply(lambda x: x.sort_values('VAF%_flt3_ext',
			ascending=False))
			
			merged_df = merged_df.sort_values(by=['VAF_getITD(%)', 'VAF_FiLT3r(%)', 'VAF%_flt3_ext'], ascending=[False, False, False])

		except Exception as e:
			print(f"❌ Error during merging: {e}")
			merged_df = pd.DataFrame()
			proceed_merge = False

#save outfile
with pd.ExcelWriter(args.output) as writer:
	merged_df.to_excel(writer, sheet_name = args.sample_id + "_common", index=False)
	df_flt3_itd_ext_in.to_excel(writer, sheet_name= args.sample_id +"_flt3_itd_ext", index=False)
	df_getITD_in.to_excel(writer, sheet_name= args.sample_id + "_getITD", index=False)
	df_flt3r_in.to_excel(writer, sheet_name= args.sample_id + "_FiLT3r", index=False)
	flt3r_df.to_excel(writer, sheet_name= args.sample_id +"_Filt3r_anno", index=False)
	varscan_df.to_excel(writer, sheet_name= args.sample_id + "_varscan", index=False)
	cov_bed.to_excel(writer, sheet_name="coverage", index=False)



