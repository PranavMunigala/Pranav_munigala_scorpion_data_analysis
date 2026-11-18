import pandas as pd, numpy as np
pd.set_option('display.width',200); pd.set_option('display.max_columns',50); pd.set_option('display.max_rows',250)
D='dataset/'
d=pd.read_csv(D+'design.csv',dtype={'uuid':str})
print("rows:",len(d),"unique uuid:",d.uuid.nunique())
print("\n--- design_set / format ---")
print(d.groupby(['design_set','intended_toxins','format']).size())
print("\n--- CLAIM: 211 designs, groups 60/30/30/59+32 ---")
print(d.intended_toxins.value_counts())
print(d.format.value_counts())
print("\n--- binder_against counts (CLAIM 17: Cn2 7, AaHII 4, AaHI 4, LqhIII 2, Ts1 0) ---")
b=d[d.binder_against.notna() & (d.binder_against.astype(str).str.strip()!='')]
print("n binders:",len(b))
print(b.binder_against.value_counts())
print("\n--- the 17 binders, KD table (compare to Table 1) ---")
toxmap={'Cn2':'Cn2','Ts1':'Ts1','AaH II':'AaHII','AaH I':'AaHI','LqhIII':'LqhIII'}
rows=[]
for _,r in b.iterrows():
    t=r.binder_against.strip(); k=toxmap.get(t,t)
    rows.append(dict(name=r['name'],uuid=r.uuid,tox=t,made_against=r.intended_toxins,fmt=r.format,len=r.length_aa,
        adaptyv=r.get(f'{k}_Adaptyv_KD_nM'), twist=r.get(f'{k}_Twist_KD_nM') if f'{k}_Twist_KD_nM' in d.columns else None))
bt=pd.DataFrame(rows).sort_values(['tox','twist','adaptyv'])
print(bt.to_string(index=False))
print("\n--- lab KD ratio (where both) ---")
bt2=bt.dropna(subset=['adaptyv','twist']).copy()
bt2['ratio']=(bt2.adaptyv/bt2.twist)
print(bt2[['name','tox','adaptyv','twist','ratio']].to_string(index=False))
print("fold-range overall: %.2f to %.2f"%(bt2.ratio.min(),bt2.ratio.max()))
