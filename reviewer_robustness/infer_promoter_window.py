import gzip,re,bisect,pandas as pd,numpy as np
from collections import defaultdict
from pathlib import Path
base=Path('/mnt/data')
tx=defaultdict(list)
with gzip.open(base/'mm10.refGene.gtf.gz','rt') as f:
 for line in f:
  if line.startswith('#'): continue
  p=line.rstrip().split('\t')
  if len(p)<9 or p[2]!='transcript': continue
  ch=p[0]; s=int(p[3])-1; e=int(p[4]); strand=p[6]
  tx[ch].append((s,e,strand))

def merge(iv):
 if not iv:return []
 iv=sorted(iv);out=[list(iv[0])]
 for s,e in iv[1:]:
  if s<=out[-1][1]:out[-1][1]=max(out[-1][1],e)
  else:out.append([s,e])
 return [tuple(x) for x in out]
body={ch:merge([(s,e) for s,e,_ in arr]) for ch,arr in tx.items()}
def make_prom(up,down):
 d={}
 for ch,arr in tx.items():
  iv=[]
  for s,e,strand in arr:
   t=s if strand=='+' else e-1
   if strand=='+': ps=max(0,t-up);pe=t+down+1
   else: ps=max(0,t-down);pe=t+up+1
   iv.append((ps,pe))
  d[ch]=merge(iv)
 return d
def inside(iv,p):
 starts=[x[0] for x in iv];i=bisect.bisect_right(starts,p)-1;return i>=0 and p<iv[i][1]
D=pd.concat([pd.read_csv(base/'attached_blood_only_annotated.csv'),pd.read_csv(base/'attached_liver_only_annotated.csv'),pd.read_csv(base/'attached_shared_annotated.csv')],ignore_index=True).drop_duplicates(['chr','pos'])
obs=D.gene_region.astype(str).to_numpy()
for up,down in [(500,500),(1000,500),(1500,500),(2000,500),(2000,1000),(2000,2000),(1000,1000),(5000,500)]:
 pr=make_prom(up,down); starts={ch:[x[0] for x in iv] for ch,iv in pr.items()}; bstarts={ch:[x[0] for x in iv] for ch,iv in body.items()}
 pred=[]
 for ch,pos in zip(D.chr.astype(str),D.pos.astype(int)):
  def hit(di,st):
   iv=di.get(ch,[]);ss=st.get(ch,[]);i=bisect.bisect_right(ss,pos)-1;return i>=0 and pos<iv[i][1]
  pred.append('promoter' if hit(pr,starts) else ('gene_body' if hit(body,bstarts) else 'intergenic'))
 pred=np.array(pred)
 acc=(pred==obs).mean(); prom_recall=((pred=='promoter')&(obs=='promoter')).sum()/max(1,(obs=='promoter').sum()); prom_prec=((pred=='promoter')&(obs=='promoter')).sum()/max(1,(pred=='promoter').sum())
 print(up,down,'acc',acc,'prom recall',prom_recall,'prec',prom_prec,'predprom',sum(pred=='promoter'),'obsprom',sum(obs=='promoter'))
