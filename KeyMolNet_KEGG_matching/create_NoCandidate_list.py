import pandas as pd
import sys

all_df = pd.read_table("../matching_output/KeyMolNet_icd10_noItems.txt")
## add ICDclass to all_df
class_df = pd.read_csv("../disease_raw/Disease ID_KeyMolnet.csv")
print(class_df)
## add KMN_ID to class_df
class_df["KMN_ID"] = [ f"KMN{x}" for x in class_df["No"]]
class_df = class_df.set_index("KMN_ID")

all_df = all_df.set_index("KMN_ID")
print(all_df)

all_df = pd.concat([all_df, class_df.loc[:, ["ICD___", "Disease_english"]]], axis=1)
all_df = all_df.rename(columns={"ICD___":"ICD_class"})
print(all_df)
## save all_df
# all_df.to_csv("../matching_output/KMN_wICDcode_0714.txt", sep="\t")

##
candidate_df = pd.read_table("KMN_candidate_list.txt")
print(candidate_df)
candidate_list = list(candidate_df["KMN_candidate_IDs"])

all_df = all_df.reset_index()
print(all_df)

NoCandidate_df = all_df[~all_df["index"].isin(candidate_list)]
NoCandidate_df = NoCandidate_df.set_index("index").drop(columns=["Count", "Name"]).reindex(columns=["Disease_english", "ICD_class", "ICD10"])
print(NoCandidate_df)
## save NoCandidate_df
# NoCandidate_df.to_csv("KMN_462diseases_0714.csv")

## add all ICD10 codes to all_df
all_df = all_df.set_index("index")
print(all_df)
final_0714_df = pd.read_csv("../KeyMolNet_KEGG_matching/KMN_with_ICD10code_0723.csv", index_col=0)
print(final_0714_df)

for kmn in final_0714_df.index:
    print(kmn)
    ## updated icd10 code in all_df
    all_df.at[kmn, "ICD10"] = final_0714_df.at[kmn, "ICD10"]

print(all_df)
all_df.to_csv("../matching_output/KMN_wICDcode_0723.txt", sep="\t")
