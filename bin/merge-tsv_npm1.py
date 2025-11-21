#!/usr/bin/env python3
import pandas as pd
import os
import sys

output = sys.argv[1]
varscan_file = sys.argv[2]
coverage_file = sys.argv[3]

writer = pd.ExcelWriter(output, engine='openpyxl')

if os.path.getsize(varscan_file) != 0:
    df_varscan = pd.read_csv(varscan_file, sep=',')
    if not df_varscan.empty:
        df_varscan.to_excel(writer, sheet_name='varscan', index=False)
    else:
        pd.DataFrame(columns=df_varscan.columns).to_excel(writer, sheet_name='varscan', index=False)
else:
    print("Skipping varscan (empty file).")

if os.path.getsize(coverage_file) != 0:
    df_coverage = pd.read_csv(coverage_file, sep='\t')
    if not df_coverage.empty:
        df_coverage.to_excel(writer, sheet_name='coverage', index=False, header=False)
    else:
        pd.DataFrame(columns=df_coverage.columns).to_excel(writer, sheet_name='coverage', index=False, header=False)
else:
    print("Skipping coverage (empty file).")

writer.close()

