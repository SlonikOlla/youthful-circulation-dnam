import pandas as pd, numpy as np
from pathlib import Path
base=Path('/mnt/data')
b=pd.read_csv(base/'attached_blood_only_annotated.csv')
l=pd.read_csv(base/'attached_liver_only_annotated.csv')
b['B']=1;b['L']=0;l['B']=0;l['L']=1
u=pd.concat([b,l],ignore_index=True).drop_duplicates(['chr','pos']).copy()
canon=[f'chr{i}' for i in range(1,20)]+['chrX','chrY']
u=u[u.chr.isin(canon)].copy()
genic=u.gene_region.isin(['promoter','gene_body']) & u.mapped_gene.notna() & (u.mapped_gene.astype(str)!='')
genes=sorted(u.loc[genic,'mapped_gene'].astype(str).unique()); gid={g:i for i,g in enumerate(genes)}
u['gid']=-1
u.loc[genic,'gid']=u.loc[genic,'mapped_gene'].astype(str).map(gid).astype(int)
bg=set(u.loc[(u.B==1)&(u.gid>=0),'gid']); lg=set(u.loc[(u.L==1)&(u.gid>=0),'gid'])
obs=len(bg&lg)
print('observed canonical genic overlap',obs,'B genes',len(bg),'L genes',len(lg),'total rows',len(u))
blocks=[]
for ch,d in u.sort_values(['chr','pos']).groupby('chr',sort=False):
 blocks.append((ch,d.B.to_numpy(np.int8),d.L.to_numpy(np.int8),d.gid.to_numpy(np.int32)))
rng=np.random.default_rng(20260819)
N=5000
null=np.empty(N,dtype=np.int32)
for it in range(N):
 bgenes=set(); lgenes=set()
 for ch,B,L,G in blocks:
  n=len(G)
  sb=int(rng.integers(0,n)); sl=int(rng.integers(0,n))
  br=np.roll(B,sb).astype(bool); lr=np.roll(L,sl).astype(bool)
  gb=G[br & (G>=0)]; gl=G[lr & (G>=0)]
  if len(gb): bgenes.update(np.unique(gb).tolist())
  if len(gl): lgenes.update(np.unique(gl).tolist())
 null[it]=len(bgenes & lgenes)
mean=null.mean(); lo,hi=np.quantile(null,[.025,.975]); p=(1+(null>=obs).sum())/(N+1)
print('null mean',mean,'95%',lo,hi,'fold',obs/mean,'p',p,'sd',null.std(ddof=1))
pd.DataFrame({'iteration':np.arange(1,N+1),'common_genes':null}).to_csv(base/'genic_promoter_genebody_null_distribution.csv',index=False)
pd.DataFrame([{'observed_common_genes':obs,'blood_genes':len(bg),'liver_genes':len(lg),'null_mean':mean,'null_sd':null.std(ddof=1),'null_2.5':lo,'null_97.5':hi,'fold_vs_null':obs/mean,'empirical_p':p,'iterations':N,'seed':20260819,'null_definition':'independent chromosome-wise circular shifts of blood/liver labels across union of tissue-specific reversal loci; canonical chromosomes; promoter+gene_body only'}]).to_csv(base/'genic_promoter_genebody_circular_null_summary.csv',index=False)
