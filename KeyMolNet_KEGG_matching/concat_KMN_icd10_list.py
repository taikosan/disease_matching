import pandas as pd


df1 = pd.read_csv("KMN_candidate_list_filled.csv")
df1 = df1.rename(columns={"KMN_candidate_IDs":"index", "KMN_candidate_Names":"Disease_english", "ICDclass":"ICD_class", "KMN_ICD10":"ICD10"})


df2 = pd.read_csv("KMN_462diseases_0714_final.csv")

print(df1)
print(df2)

all_df = pd.concat([df1, df2], axis=0)
print(all_df)

all_df.to_csv("KMN_with_ICD10code_0723.csv", index=False)
