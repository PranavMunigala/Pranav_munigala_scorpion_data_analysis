import pandas as pd, numpy as np
pd.set_option('display.width',250);pd.set_option('display.max_rows',300)
d=pd.read_csv('dataset/design.csv',dtype={'uuid':str})
s=pd.read_csv('dataset/in_silico_scores.csv',dtype={'uuid':str})
o=pd.read_csv('dataset/supporting_tables/origin.csv',dtype={'uuid':str,'parent_uuid':str,'domainA_uuid':str,'domainB_uuid':str})
P=pd.read_csv('analysis/pairs_final.csv',dtype={'uuid':str})
TOXK={'Cn2':'Cn2','Ts1':'Ts1','AaH II':'AaHII','AaH I':'AaHI','LqhIII':'LqhIII'}
P['k']=P.tox.map(TOXK)
M=P.merge(s,on='uuid',how='left',suffixes=('','_s'))
M['final']=[r[f'{r.k}_final_score'] for _,r in M.iterrows()]
M['ipsae_mean']=[np.nanmean([r[f'{r.k}_ipSAE_ESMFold2'],r[f'{r.k}_ipSAE_ESMFold2Fast'],r[f'{r.k}_ipSAE_ProtenixV2']]) for _,r in M.iterrows()]
M['scdockq_mean']=[np.nanmean([r[f'{r.k}_scDockQ_ESMFold2'],r[f'{r.k}_scDockQ_ESMFold2Fast'],r[f'{r.k}_scDockQ_ProtenixV2']]) for _,r in M.iterrows()]

def auc(y,x):
    y=np.asarray(y,bool); x=np.asarray(x,float)
    m=~np.isnan(x); y,x=y[m],x[m]
    if y.sum()==0 or (~y).sum()==0: return np.nan,y.sum(),(~y).sum()
    r=pd.Series(x).rank().values
    a=(r[y].sum()-y.sum()*(y.sum()+1)/2)/(y.sum()*(~y).sum())
    return a,y.sum(),(~y).sum()

print("### IN-SILICO SCORE vs OUTCOME (AUROC; 0.5 = no discrimination) ###")
print("\n-- ALL 665 intended pairs pooled --")
for col in ['final','ipsae_mean','scdockq_mean']:
    a,np_,nn=auc(M.labelled,M[col]); print(f"  {col:13s} AUC={a:.3f}  (binders={np_}, non={nn})")

print("\n-- WITHIN the single-toxin sets only (sets 1-5, the paper's cleanest comparison) --")
S=M[M.set<=5]
for col in ['final','ipsae_mean','scdockq_mean']:
    a,np_,nn=auc(S.labelled,S[col]); print(f"  {col:13s} AUC={a:.3f}  (binders={np_}, non={nn})")

print("\n-- WITHIN EACH single-toxin set (n=12 each) --")
for st,t in zip(range(1,6),['Cn2','Ts1','AaH II','AaH I','LqhIII']):
    sub=M[M.set==st]
    a,np_,nn=auc(sub.labelled,sub['final'])
    ai,_,_=auc(sub.labelled,sub['ipsae_mean'])
    print(f"  set{st} {t:8s} n={len(sub):2d} binders={np_}  AUC(final)={a if np.isnan(a) else round(a,3)}  AUC(ipSAE)={ai if np.isnan(ai) else round(ai,3)}")

print("\n### MEDIAN final score BY GROUP vs observed success rate (the paper's 'tracks across sets' claim) ###")
M['grp']=np.where(M.set<=5,'1 single-toxin',np.where(M.set==6,'2 double beta',np.where(M.set==7,'3 triple alpha','4 all-five')))
g=M.groupby('grp').agg(pairs=('uuid','size'),median_final=('final','median'),median_ipsae=('ipsae_mean','median'))
bd=d.assign(isb=d.binder_against.notna()&(d.binder_against.astype(str).str.strip()!=''))
bd['grp']=np.where(bd.design_set<=5,'1 single-toxin',np.where(bd.design_set==6,'2 double beta',np.where(bd.design_set==7,'3 triple alpha','4 all-five')))
g2=bd.groupby('grp').agg(designs=('uuid','size'),binders=('isb','sum'))
g2['success_%']=(100*g2.binders/g2.designs).round(1)
print(g.join(g2).to_string())
print("\nSpearman(median_final, success_%) across the 4 groups:",
      round(pd.Series(pd.Series(g.median_final.values).rank().values).corr(pd.Series(pd.Series(g2["success_%"].values).rank().values)),3))

print("\n### RANK-IN-SET of the binders (Claude's own ranking: 1 = its top pick) ###")
R=d.merge(o[['uuid','rank_in_set','design_tool','backbone_family','backbone_made_against','sequence_design_method']],on='uuid',how='left')
R['isb']=R.binder_against.notna()&(R.binder_against.astype(str).str.strip()!='')
for st,t in zip(range(1,6),['Cn2','Ts1','AaH II','AaH I','LqhIII']):
    sub=R[R.design_set==st].sort_values('rank_in_set')
    hits=sub[sub.isb].rank_in_set.tolist()
    print(f"  set{st} {t:8s} binder ranks among 12: {hits}")
print("\n  all-five + multi sets binder ranks:")
for st in [6,7,8,9,10]:
    sub=R[R.design_set==st]; print(f"    set{st}: n={len(sub)} binder ranks={sub[sub.isb].rank_in_set.tolist()}")
print("\n### DESIGN TOOL attribution ###")
tt=R.groupby('design_tool').agg(designs=('uuid','size'),binders=('isb','sum'))
tt['rate_%']=(100*tt.binders/tt.designs).round(1)
print(tt.sort_values('binders',ascending=False).to_string())
print("\n### sequence design method ###")
print(R.groupby('sequence_design_method').agg(designs=('uuid','size'),binders=('isb','sum')).to_string())
M.to_csv('analysis/pairs_scored.csv',index=False)
