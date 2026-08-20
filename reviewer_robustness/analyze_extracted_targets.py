#!/usr/bin/env python3
import argparse, pandas as pd, numpy as np
from scipy.stats import pearsonr, spearmanr

p=argparse.ArgumentParser()
p.add_argument('input')
p.add_argument('--out',default='replicate_aware_DMC_summary.csv')
a=p.parse_args()
df=pd.read_csv(a.input)
# This utility summarizes replicate-aware locus-level statistics exported from the analysis workflow.
# Expected columns include tissue, comparison/effect fields, and locus-level P/FDR values.
summary=[]
for tissue,g in df.groupby('tissue') if 'tissue' in df.columns else [('all',df)]:
    row={'tissue':tissue,'n_rows':len(g)}
    for col in g.columns:
        lc=col.lower()
        if 'fdr' in lc or 'qvalue' in lc or lc in ('padj','adj_p'):
            x=pd.to_numeric(g[col],errors='coerce')
            row[f'{col}_lt_0.05']=int((x<0.05).sum())
    summary.append(row)
pd.DataFrame(summary).to_csv(a.out,index=False)
print(pd.DataFrame(summary).to_string(index=False))
