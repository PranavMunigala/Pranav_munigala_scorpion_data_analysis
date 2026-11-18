import pandas as pd, numpy as np
pd.set_option('display.width',250);pd.set_option('display.max_rows',400);pd.set_option('display.max_columns',40)
s=pd.read_csv('dataset/sensorgrams/INDEX.csv',dtype={'uuid_or_control':str})
print("assays:",s.assay.unique()); print("labs:",s.lab.unique())
main=s[s.assay.str.contains('main',case=False,na=False)]
print("\nmain-screen rows:",len(main),"| results:",main.result.value_counts().to_dict())

print("\n=== per-replicate check on the two anomalies ===")
for u in ['9a8cdc7a','cf8652bc']:
    sub=main[(main.uuid_or_control==u)]
    print(f"\n--- {u} ---")
    print(sub[['lab','toxin','replicate','result','assay','format']].to_string(index=False))

print("\n=== CRITERION 1 re-derived properly: hit in EVERY replicate, both labs ===")
d=pd.read_csv('dataset/design.csv',dtype={'uuid':str})
TOXK={'Cn2':'Cn2','Ts1':'Ts1','AaH II':'AaHII','AaH I':'AaHI','LqhIII':'LqhIII'}
out=[]
for _,r in d.iterrows():
    for t in [x.strip() for x in str(r.intended_toxins).split(';')]:
        k=TOXK[t]
        reps=main[(main.uuid_or_control==r.uuid)&(main.toxin==t)]
        labs=reps.lab.unique()
        everyrep = len(reps)>0 and all(str(x)=='hit' for x in reps.result)
        akd=r.get(f'{k}_Adaptyv_KD_nM'); tkd=r.get(f'{k}_Twist_KD_nM') if f'{k}_Twist_KD_nM' in d else np.nan
        rmax=r.get(f'{k}_Twist_Rmax_ratio') if f'{k}_Twist_Rmax_ratio' in d else np.nan
        kds=[x for x in [akd,tkd] if pd.notna(x)]
        c1 = everyrep and (len(labs)>=2 if t!='LqhIII' else True)
        c2 = bool(kds) and all(x<=1000 for x in kds)
        c3 = True if t=='LqhIII' else (pd.notna(rmax) and rmax<=1.0)
        c4 = (str(r.CLEC12A_Twist_result_run_1)!='hit') and (str(r.CLEC12A_Twist_result_run_2)!='hit')
        out.append(dict(uuid=r.uuid,name=r['name'],tox=t,nrep=len(reps),labs=len(labs),
            c1=c1,c2=c2,c3=c3,c4=c4,allp=c1 and c2 and c3 and c4,
            labelled=(str(r.binder_against).strip()==t)))
O=pd.DataFrame(out)
print("my passes:",O.allp.sum()," labelled:",O.labelled.sum())
dis=O[O.allp!=O.labelled]
print("disagreements:",len(dis))
print(dis.to_string(index=False) if len(dis) else "  NONE — exact reproduction of all 17 binder labels")
