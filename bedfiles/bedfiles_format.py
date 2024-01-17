#!/usr/bin/env python3

import sys, csv

infile = sys.argv[1]
outfile = sys.argv[2]

output_file = open (outfile,'w')
with open (infile, 'r') as input_file:
	tsv_handle = csv.reader(input_file, delimiter="\t")
	for lines in tsv_handle:
		#gene_name = lines[3].split(';')[2]
		#gene_name=str(lines[0]) + '_' +''.join(lines[1]) + '_' +''.join(lines[2])
		chrom = str(lines[0]).split(':')[0]
		start = str(lines[0]).split(':')[1].split('-')[0]
		end = str(lines[0]).split(':')[1].split('-')[1]
		region_name = lines[1]
		print (chrom, start, end, region_name, file=output_file, sep='\t')
		#print (lines[0], lines[1], lines[2], gene_name, file=output_file, sep='\t')
		#print (lines[0], lines[1], lines[2], lines[3], lines[4], file=outfile, sep='\t')

output_file.close()
