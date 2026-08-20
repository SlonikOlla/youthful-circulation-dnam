import pandas as pd, numpy as np, glob, os
from scipy.stats import ttest_ind
from pathlib import Path
base=Path('/mnt/data')
mapdf=pd.read_csv(base/'sample_map.tsv',sep='\t',names=['gsm','tissue','group'])

def bh_full(p,m_total=1008133):
    p=np.asarray(p,float); order=np.argsort(p); ps=p[order]
    q=ps*m_total/np.arange(1,len(p)+1)
    q=np.minimum.accumulate(q[::-1])[::-1]; q=np.clip(q,0,1)
    out=np.empty_like(q); out[order]=q; return out
rows=[]; details=[]
for tissue in ['blood','liver']:
    if tissue=='blood':
        targ=pd.concat([pd.read_csv(base/'attached_blood_only_annotated.csv'),pd.read_csv(base/'attached_shared_annotated.csv')],ignore_index=True)
    else:
        targ=pd.concat([pd.read_csv(base/'attached_liver_only_annotated.csv'),pd.read_csv(base/'attached_shared_annotated.csv')],ignore_index=True)
    targ=targ.drop_duplicates(['chr','pos']).sort_values(['chr','pos']).reset_index(drop=True)
    keys=list(zip(targ.chr.astype(str),targ.pos.astype(int))); index={k:i for i,k in enumerate(keys)}
    mats={g:[] for g in ['YoungIso','OldIso','OldHet']}
    covmats={g:[] for g in mats}
    for _,r in mapdf[mapdf.tissue==tissue].iterrows():
        fn=glob.glob(str(base/'extracted_targets'/f"{r.gsm}_{tissue}_{r.group}.tsv"))[0]
        d=pd.read_csv(fn,sep='\t',names=['chr','pos','meth','unmeth'])
        a=np.empty(len(keys),dtype=float); c=np.empty(len(keys),dtype=float)
        for ch,pos,m,u in d.itertuples(index=False):
            i=index[(str(ch),int(pos))]; tot=m+u; a[i]=m/tot; c[i]=tot
        mats[r.group].append(a); covmats[r.group].append(c)
    Y,O,H=[np.vstack(mats[g]) for g in ['YoungIso','OldIso','OldHet']]
    pa=ttest_ind(O,Y,axis=0,equal_var=False).pvalue; ph=ttest_ind(H,O,axis=0,equal_var=False).pvalue
    qa=bh_full(pa); qh=bh_full(ph)
    age=O.mean(0)-Y.mean(0); hpb=H.mean(0)-O.mean(0)
    rev=(np.sign(age)==-np.sign(hpb))&(np.abs(age)>=.10)&(np.abs(hpb)>=.10)
    for lab,mask in [('effect10',rev),('both_p05',rev&(pa<.05)&(ph<.05)),('both_q10',rev&(qa<.10)&(qh<.10)),('both_q05',rev&(qa<.05)&(qh<.05))]:
        rows.append([tissue,lab,int(mask.sum()),mask.sum()/rev.sum()])
    out=targ[['chr','pos']].copy(); out['age_delta']=age; out['hpb_delta']=hpb; out['age_p']=pa; out['hpb_p']=ph; out['age_q_conservative']=qa; out['hpb_q_conservative']=qh; out['effect10_reversal']=rev; out['both_q05']=rev&(qa<.05)&(qh<.05)
    out.to_csv(base/f'{tissue}_replicate_aware_effect10_stats.csv',index=False)
pd.DataFrame(rows,columns=['tissue','criterion','n','fraction_of_effect10']).to_csv(base/'replicate_aware_DMC_summary.csv',index=False)
print(pd.DataFrame(rows,columns=['tissue','criterion','n','fraction_of_effect10']).to_string(index=False))
