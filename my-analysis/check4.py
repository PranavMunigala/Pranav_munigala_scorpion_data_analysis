import pandas as pd, numpy as np
pd.set_option('display.width',250);pd.set_option('display.max_rows',300)
d=pd.read_csv('dataset/design.csv',dtype={'uuid':str})
print("CLEC12A run1 values:",d.CLEC12A_Twist_result_run_1.value_counts(dropna=False).to_dict())
print("CLEC12A run2 values:",d.CLEC12A_Twist_result_run_2.value_counts(dropna=False).to_dict())
b=d[d.binder_against.notna()&(d.binder_against.astype(str).str.strip()!='')]
print("\nCLEC12A calls among the 17 LABELLED binders:")
print(b[['name','uuid','CLEC12A_Twist_result_run_1','CLEC12A_Twist_result_run_2']].to_string(index=False))
print("\nAll designs whose CLEC12A has 'near background' in either run:")
nb=d[(d.CLEC12A_Twist_result_run_1=='near background')|(d.CLEC12A_Twist_result_run_2=='near background')]
print(nb[['uuid','name','binder_against','CLEC12A_Twist_result_run_1','CLEC12A_Twist_result_run_2']].to_string(index=False))
print("\n--- KD_note for the two anomalies & binders ---")
for u in ['9a8cdc7a','cf8652bc']:
    r=d[d.uuid==u].iloc[0]
    print(u,"| binder_against=",repr(r.binder_against),"| KD_note=",repr(r.KD_note))
