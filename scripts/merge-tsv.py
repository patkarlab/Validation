#!/usr/bin/env python3
import pandas as pd
import os
import sys

args = sys.argv
output = args[1]  

tsvfilenames = ['varscan', 'filt3r']
filenameindex = 0

writer = pd.ExcelWriter(output)  

for arguments in args[2:]:  
    if os.path.getsize(arguments) != 0:  # File is not empty
        df = pd.read_csv(arguments, sep=',')
        if not df.empty:  # File contains data
            df.to_excel(writer, sheet_name=tsvfilenames[filenameindex], index=False)
        else:  # File has only the header
            print(f"Adding sheet for {tsvfilenames[filenameindex]} with only header.")
            pd.DataFrame(columns=df.columns).to_excel(
                writer, sheet_name=tsvfilenames[filenameindex], index=False
            )
    else:  # File has no content (neither header nor data)
        print(f"Skipping {tsvfilenames[filenameindex]} as it has no content.")
    filenameindex += 1

writer.save()
print(f"Output written to {output}")
