import pandas as pd
import sys


df = pd.read_csv("icd10cm_order_2021.csv", header=None)
df.columns = ["order", "raw_code", "hatena", "Name1", "Name2"]
print(df)


raw_code_list = list(df["raw_code"])
print(len(raw_code_list))

icd10_code_list = []
for code in raw_code_list:
    if len(code) == 3:
        new_code = code
    else:
        new_code = f"{code[0:3]}.{code[3:]}"
    icd10_code_list.append(new_code)


print(len(icd10_code_list))

df["icd10_code"] = icd10_code_list

df = df.reindex(columns=["order", "raw_code", "icd10_code", "Name1", "hatena", "Name2"])
print(df)

df_new = df.loc[:, ["order", "raw_code", "icd10_code", "Name1"]].set_index("order")
print(df_new)

df_new.to_csv("ICD10_code_Name_2021.csv")
