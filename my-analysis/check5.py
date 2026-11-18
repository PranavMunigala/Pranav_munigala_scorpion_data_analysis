import pandas as pd, numpy as np, glob, os
pd.set_option('display.width',250);pd.set_option('display.max_rows',300)
d=pd.read_csv('dataset/design.csv',dtype={'uuid':str})
K=[]
for f in glob.glob('dataset/sensorgrams/*/*_kinetics.csv'):
    u=os.path.basename(f).split('_')[0]
    try: k=pd.read_csv(f)
    except Exception: continue
    k['uuid']=u; K.append(k)
K=pd.concat(K,ignore_index=True)
print("kinetics rows:",len(K),"designs/controls with folders:",K.uuid.nunique())
print("cols:",list(K.columns))

TOXK={'Cn2':'Cn2','Ts1':'Ts1','AaH II':'AaHII','AaH I':'AaHI','LqhIII':'LqhIII'}
def percall(u,tox):
    """Adaptyv per-replicate calls (prefer Adaptyv_replicate_result), + Twist call"""
    sub=K[(K.uuid==u)&(K.toxin==tox)]
    a=sub[sub.lab=='Adaptyv Bio']; t=sub[sub.lab=='Twist Bioscience']
    acalls=[]
    for _,r in a.iterrows():
        v=r.get('Adaptyv_replicate_result')
        acalls.append(str(v) if pd.notna(v) else str(r.result))
    tcalls=[str(x) for x in t.result]
    return acalls,tcalls

out=[]
for _,r in d.iterrows():
    for tox in [x.strip() for x in str(r.intended_toxins).split(';')]:
        k=TOXK[tox]
        ac,tc=percall(r.uuid,tox)
        akd=r.get(f'{k}_Adaptyv_KD_nM'); tkd=r.get(f'{k}_Twist_KD_nM') if f'{k}_Twist_KD_nM' in d.columns else np.nan
        rmax=r.get(f'{k}_Twist_Rmax_ratio') if f'{k}_Twist_Rmax_ratio' in d.columns else np.nan
        kds=[x for x in [akd,tkd] if pd.notna(x)]
        if tox=='LqhIII':
            c1 = len(ac)>0 and all(x=='hit' for x in ac)
        else:
            c1 = len(ac)>0 and all(x=='hit' for x in ac) and len(tc)>0 and all(x=='hit' for x in tc)
        c2 = bool(kds) and all(x<=1000 for x in kds)
        c3 = True if tox=='LqhIII' else (pd.notna(rmax) and rmax<=1.0)
        c4 = (str(r.CLEC12A_Twist_result_run_1)=='no hit') and (str(r.CLEC12A_Twist_result_run_2)=='no hit')
        out.append(dict(uuid=r.uuid,name=r['name'],tox=tox,set=r.design_set,
                        arep='/'.join(ac),trep='/'.join(tc),akd=akd,tkd=tkd,rmax=rmax,
                        c1=c1,c2=c2,c3=c3,c4=c4,allp=all([c1,c2,c3,c4]),
                        labelled=(str(r.binder_against).strip()==tox)))
O=pd.DataFrame(out)
print("\n### DEFINITIVE REPRODUCTION ###")
print("pairs:",len(O)," my passes:",O.allp.sum()," paper's labels:",O.labelled.sum())
dis=O[O.allp!=O.labelled]
print("disagreements:",len(dis))
if len(dis): print(dis[['uuid','name','tox','arep','trep','akd','tkd','rmax','c1','c2','c3','c4','allp','labelled']].to_string(index=False))
else: print(">>> EXACT MATCH: all 17 binder labels reproduced from the raw released data.")
O.to_csv('analysis/pairs_final.csv',index=False)

# success rates by group
print("\n### SUCCESS RATES BY DESIGN GROUP (paper: 11/60=18% single, 5/60=8% double+triple, 1/91 all-five) ###")
d['grp']=np.where(d.design_set<=5,'single-toxin (sets 1-5)',
          np.where(d.design_set==6,'double beta (set 6)',
          np.where(d.design_set==7,'triple alpha (set 7)','all-five (sets 8-10)')))
bb=d.assign(isb=d.binder_against.notna()&(d.binder_against.astype(str).str.strip()!=''))
g=bb.groupby('grp').agg(designs=('uuid','size'),binders=('isb','sum'))
g['rate_%']=(100*g.binders/g.designs).round(1)
print(g.to_string())
print("\nsingle+double+triple combined:",int(bb[bb.design_set<=7].isb.sum()),"/",len(bb[bb.design_set<=7]))
print("double+triple only:",int(bb[bb.design_set.isin([6,7])].isb.sum()),"/",len(bb[bb.design_set.isin([6,7])]))
print("per-toxin single-toxin sets:")
for s,t in zip(range(1,6),['Cn2','Ts1','AaH II','AaH I','LqhIII']):
    sub=bb[bb.design_set==s]; print(f"  set{s} {t:8s} {int(sub.isb.sum())}/{len(sub)}")
