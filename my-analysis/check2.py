import pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_columns',60); pd.set_option('display.max_rows',300)
d=pd.read_csv('dataset/design.csv',dtype={'uuid':str})
TOX=[('Cn2','Cn2'),('Ts1','Ts1'),('AaH II','AaHII'),('AaH I','AaHI'),('LqhIII','LqhIII')]
def intended(r,t): return t in [x.strip() for x in str(r.intended_toxins).split(';')]

print("="*70); print("A. HIT TALLY (main screen) per toxin / lab")
for t,k in TOX:
    sub=d[d.apply(lambda r: intended(r,t),axis=1)]
    a=sub[f'{k}_Adaptyv_result'].value_counts().to_dict()
    tw=sub[f'{k}_Twist_result'].value_counts().to_dict() if f'{k}_Twist_result' in d else {}
    print(f"{t:8s} intended={len(sub):3d}  Adaptyv={a}  Twist={tw}")

print("="*70); print("B. RE-DERIVE the 4 binder criteria from raw columns")
recon=[]
for _,r in d.iterrows():
    for t,k in TOX:
        if not intended(r,t): continue
        ar=str(r.get(f'{k}_Adaptyv_result','')); tr=str(r.get(f'{k}_Twist_result','')) if f'{k}_Twist_result' in d else 'nm'
        akd=r.get(f'{k}_Adaptyv_KD_nM'); tkd=r.get(f'{k}_Twist_KD_nM') if f'{k}_Twist_KD_nM' in d else np.nan
        rmax=r.get(f'{k}_Twist_Rmax_ratio') if f'{k}_Twist_Rmax_ratio' in d else np.nan
        c1=(ar=='hit') and (tr=='hit' or t=='LqhIII')
        kds=[x for x in [akd,tkd] if pd.notna(x)]
        c2=bool(kds) and all(x<=1000 for x in kds)
        c3=True if t=='LqhIII' else (pd.notna(rmax) and rmax<=1.0)
        c4=(str(r.CLEC12A_Twist_result_run_1)!='hit') and (str(r.CLEC12A_Twist_result_run_2)!='hit')
        recon.append(dict(name=r['name'],uuid=r.uuid,tox=t,set=r.design_set,fmt=r.format,
            a=ar,tw=tr,akd=akd,tkd=tkd,rmax=rmax,c1=c1,c2=c2,c3=c3,c4=c4,
            all4=c1 and c2 and c3 and c4, labelled=(str(r.binder_against).strip()==t)))
R=pd.DataFrame(recon)
print("intended design-toxin pairs:",len(R))
print("hit at BOTH labs (Cn2/Ts1/AaHII/AaHI only):",(R[R.tox!='LqhIII'].a.eq('hit')&R[R.tox!='LqhIII'].tw.eq('hit')).sum(),
      " | LqhIII hits at Adaptyv:",(R[R.tox=='LqhIII'].a=='hit').sum())
print("\nmy all-4-criteria passes:",R.all4.sum()," | labelled binders:",R.labelled.sum())
print("\nDISAGREEMENTS between my re-derivation and their label:")
dis=R[R.all4!=R.labelled]
print(dis[['name','uuid','tox','set','a','tw','akd','tkd','rmax','c1','c2','c3','c4','all4','labelled']].to_string(index=False) if len(dis) else "  NONE - perfect reproduction")

print("="*70); print("C. WHY the both-lab hits failed (paper: 35 pairs hit both labs, 15 passed, 20 failed; 15 of those on CLEC12A)")
bl=R[((R.a=='hit')&(R.tw=='hit'))|((R.tox=='LqhIII')&(R.a=='hit'))]
print("pairs hit at both labs (incl LqhIII-Adaptyv):",len(bl))
blx=R[(R.a=='hit')&(R.tw=='hit')]
print("  strictly both-lab (excl LqhIII):",len(blx),"-> pass",blx.all4.sum(),"fail",(~blx.all4).sum())
f=blx[~blx.all4]
print("  failure reasons among those:")
print("    failed c4 (CLEC12A):",(~f.c4).sum())
print("    failed c2 (KD>1uM) only:",((~f.c2)&f.c4).sum())
print("    failed c3 (Rmax) only:",((~f.c3)&f.c4&f.c2).sum())
print("    failed c1 (not every replicate - not in this table):",((~f.c1)).sum())
print(f[['name','uuid','tox','akd','tkd','rmax','c2','c3','c4']].to_string(index=False))
R.to_csv('analysis/pairs_recon.csv',index=False)
