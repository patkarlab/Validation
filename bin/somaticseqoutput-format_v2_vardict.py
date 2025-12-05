#!/usr/bin/env python3

import pandas as pd
import sys

args = sys.argv

filename = args[1]
outfile = args[2]

df = pd.read_csv(filename)
x = df['Otherinfo']
discarded_column=df.columns.get_loc('Otherinfo')
data = dict()
somatic_cols=['Chr','Start','End','Ref','Alt','FILTER','SOMATIC_FLAG','REF_COUNT','ALT_COUNT','VAF%','Func.refGene','Gene.refGene','ExonicFunc.refGene','AAChange.refGene','cosmic84','PopFreqMax']
data.setdefault('FILTER', [])
data.setdefault('VAF%', [])
data.setdefault('REF_COUNT', [])
data.setdefault('ALT_COUNT', [])
for row in x:
	#print (row)	
	rowitems=row.split('\t')
	data['FILTER'].append(rowitems[9])
	#print (rowitems[9])
	info=rowitems[10].split(';')
	formatval=rowitems[-1].split(':')
	readcounts=formatval[1]
	data['REF_COUNT'].append(formatval[3].split(',')[0])
	data['ALT_COUNT'].append(formatval[3].split(',')[1])
	vaf_calc=(float(formatval[3].split(',')[1]) / (float(formatval[3].split(',')[1]) + float(formatval[3].split(',')[0]) )) * 100
	data['VAF%'].append(round(vaf_calc, 2))

df1=df.iloc[:,:5]
#print (data.keys())
df2=pd.DataFrame(data, columns=data.keys())
df3=df.iloc[:,5:discarded_column]

horizontal_stack = pd.concat([df1, df2, df3], axis=1)
horizontal_stack['cosmic84']=horizontal_stack['cosmic84'].astype(str).replace(',' , ';')
horizontal_stack['AAChange.refGene']=horizontal_stack['AAChange.refGene'].str.replace(',' , ';')
horizontal_stack.replace(to_replace='.', value='-1', inplace=True)
horizontal_stack=horizontal_stack.reindex(columns = somatic_cols)
horizontal_stack.to_csv(outfile, index=False)