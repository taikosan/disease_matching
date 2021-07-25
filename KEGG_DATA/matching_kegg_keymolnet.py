import pandas as pd
import sys
import pprint

def get_KEGGdf_elements(kegg_df, kegg_entries_list, kegg_colname):
    target_element_list = []
    for kegg_entry in kegg_entries_list:
        target_list = txt_to_list(kegg_df.at[kegg_entry, kegg_colname])
        for elem in target_list:
            target_element_list.append(elem)

    ## target_elements list to txt
    target_elements = list_to_txt(target_element_list)

    return target_elements

def txt_to_list(txt):
    if txt == "NO_DATA":
        _list = []
    else:
        _list = txt.split(";")
    ## remove NA
    _list = [x for x in _list if not x == "NA"]
    _list = [x for x in _list if not x == "nan"]
    return _list

def list_to_txt(_list):
    if len(_list) == 0:
        txt = "NO_DATA"
    else:
        txt = ";".join(_list)
    return txt

def main(kegg_f, matching_f, KMN_ICD10_f, KMN_all_f):
    ## read KMN_all file
    KMN_all_df = pd.read_table(KMN_all_f, index_col=0)
    print(KMN_all_df)
    
    ## read kegg file
    kegg_df = pd.read_table(kegg_f, index_col=0)
    ## add kegg df to "keymolnet_id" column
    kegg_df["keymolnet_id"] = "NO_DATA"
    kegg_df["keymolnet_disease"] = "NO_DATA"
    print(kegg_df)
    
    ## create kegg_entry_name_dict
    kegg_entry_name_dict = dict(zip(list(kegg_df.index), kegg_df["Name"]))
    kegg_entry_mesh_dict = dict(zip(list(kegg_df.index), kegg_df["MeSH"]))

    ## read id matching file
    matching_df = pd.read_csv(matching_f, index_col="KEGG_ID")
    print(matching_df)
    
    ## check id matching 
    for idx in list(matching_df.index):
        keymolnet_id = matching_df.at[idx, "KeyMolnet_ID"]
        keymolnet_disease = matching_df.at[idx, "Disease_name"]
        kegg_df.at[idx, "keymolnet_id"] = f"KMN{keymolnet_id}"
        kegg_df.at[idx, "keymolnet_disease"] = keymolnet_disease
    print(kegg_df)

    #### merge KEGG and KML using ICD10 codes.
    KMN_ICD10_df = pd.read_csv(KMN_ICD10_f)
    print(KMN_ICD10_df)
    KMN_name_dict = dict(zip(KMN_ICD10_df["index"], KMN_ICD10_df["Disease_english"]))
    
    ## check KMN_candidate and already-having KEGGID KMNs
    print("##################")
    first_KMN_list = list(kegg_df[~kegg_df["keymolnet_id"].isin(["NO_DATA"])]["keymolnet_id"])
    print(len(first_KMN_list))
    first_KMN_df = kegg_df[~kegg_df["keymolnet_id"].isin(["NO_DATA"])]
    print(first_KMN_df)
    
    ### check whether ICD10 codes of KEGG and KMN are matched or not ###
    name_matched_KEGG_KMN_df = first_KMN_df.loc[:, ["Name", "ICD-10", "keymolnet_id", "keymolnet_disease"]]
    print(name_matched_KEGG_KMN_df)
    name_matched_KMN_ICD10code_df = KMN_ICD10_df.set_index("index").loc[name_matched_KEGG_KMN_df["keymolnet_id"], ["ICD10"]].reset_index().rename(columns={"ICD10":"KMN_ICD10", "index":"keymolnet_id"})
    print(name_matched_KMN_ICD10code_df)
    ## merge KEGG and KMN icd10 codes
    name_matched_KEGG_KMN_df = pd.merge(name_matched_KEGG_KMN_df, name_matched_KMN_ICD10code_df, on="keymolnet_id")
    print(name_matched_KEGG_KMN_df)
    ## check ICD10 codes are same or not
    name_matched_KEGG_KMN_df["ICD10_check"] = (name_matched_KEGG_KMN_df["ICD-10"]==name_matched_KEGG_KMN_df["KMN_ICD10"])
    print(name_matched_KEGG_KMN_df)
    print(name_matched_KEGG_KMN_df["ICD10_check"].sum())
    name_matched_KEGG_KMN_df.to_csv("NameMatched_KEGG_KMN_ICD10check.txt", sep="\t", index=False)
    ### check ICD10 codes are same or not ... END
    
    Name_matched_KMN_KEGG_dict = dict(zip(first_KMN_df["keymolnet_id"], first_KMN_df.index)) 
    print("##################")
    
    noName_match_KMN_df = KMN_ICD10_df[~KMN_ICD10_df["index"].isin(first_KMN_list)]
    print(noName_match_KMN_df)

    second_KMN_list = list(noName_match_KMN_df["index"])
    print(len(second_KMN_list))
    print("##################")
    total_KMN_list = list(set(first_KMN_list + second_KMN_list))
    print(len(total_KMN_list))
    ### check END...
    
    ## get withoutKMN_id kegg_df
    woKMN_kegg_df = kegg_df[kegg_df["keymolnet_id"].isin(["NO_DATA"]) & ~(kegg_df["ICD-10"].isin(["NO_DATA"]))]
    print(woKMN_kegg_df)
    print(woKMN_kegg_df['ICD-10'])
    KMN_ICD10_df = KMN_ICD10_df.set_index("index")
    print(KMN_ICD10_df['ICD10'])
    
    ## create KMN_icd10 dict
    KMN_ICD10_dict = dict(zip(KMN_ICD10_df.index, KMN_ICD10_df["ICD10"]))
    print(len(KMN_ICD10_dict.keys()))
    
    ## get kegg-icd10 tidy list
    KEGG_pair_list = []
    for idx, items in zip(woKMN_kegg_df.index, woKMN_kegg_df['ICD-10']):
        item_list = txt_to_list(items)
        for item in item_list:
            KEGG_pair_list.append((idx, item))
    print(len(KEGG_pair_list))    

    KEGG_tidy_df = pd.DataFrame(KEGG_pair_list)
    KEGG_tidy_df =KEGG_tidy_df.rename(columns={0:"KEGG", 1:"ICD10"})
    print(KEGG_tidy_df)

    
    ## get KMN-icd10 tidy list
    KMN_pair_list = []
    for idx, items in zip(KMN_ICD10_df.index, KMN_ICD10_df['ICD10']):
        item_list = txt_to_list(items)
        for item in item_list:
            KMN_pair_list.append((idx, item))
    print(len(KMN_pair_list))
    KMN_tidy_df = pd.DataFrame(KMN_pair_list).rename(columns={0:"KMN", 1:"ICD10"})
    print(KMN_tidy_df)
    
    ### merged df
    KMN_KEGG_merge_df = pd.merge(KEGG_tidy_df, KMN_tidy_df, on="ICD10", how="inner", indicator=True)
    print(KMN_KEGG_merge_df)
    

    ## print keymolnet == True
    print(kegg_df[~kegg_df["keymolnet_id"].isin(["NO_DATA"])])
    ## print keymolnet == True and DrugBank == True
    print(kegg_df[~kegg_df["keymolnet_id"].isin(["NO_DATA"]) & ~kegg_df["DrugBank"].isin(["NO_DATA"])])
    ## save kegg_df -- with perfect match of KMN disease name 
    # kegg_df.to_csv("KEGG_DISEASEplus_OMIM_HMDB_KMN_210623.tsv", sep="\t")

    ## add ICD10 based KMN disease matching.
    ## 2. KEGG --> KMN
    KEGG_KMN_dict = {}
    for kegg in list(kegg_df.index):
        kmn_list = sorted(list(set(KMN_KEGG_merge_df[KMN_KEGG_merge_df["KEGG"].isin([kegg])]["KMN"])))
        kmn_name_list = [KMN_name_dict[kmn] for kmn in kmn_list]
        ## update kegg_df
        if kegg_df.at[kegg, "keymolnet_id"] == "NO_DATA":
            kegg_df.at[kegg, "keymolnet_id"] = list_to_txt(kmn_list)
            kegg_df.at[kegg, "keymolnet_disease"] = list_to_txt(kmn_name_list)
    
    print(kegg_df)
    # kegg_df.to_csv("KEGG_DISEASEplus_OMIM_HMDB_KMN_210723.tsv", sep="\t")

    ### create KMN_extention_wKEGG.tsv
    ## 1. KMN --> KEGG
    KMN_KEGG_dict = {}
    for kmn in second_KMN_list:
        kegg_list = sorted(list(set(KMN_KEGG_merge_df[KMN_KEGG_merge_df["KMN"].isin([kmn])]["KEGG"])))
        KMN_KEGG_dict[kmn] = kegg_list

    KMN_extention_df = pd.DataFrame(index=KMN_all_df.index, columns=["Keymolnet_disease", "ICD10", "KEGG_ENTRY", "KEGG_Name"])
    print(KMN_extention_df)

    for kmn in list(KMN_all_df.index):
        if kmn in KMN_KEGG_dict.keys():
            kegg_entries = KMN_KEGG_dict[kmn]
            kegg_names = [kegg_entry_name_dict[kegg] for kegg in kegg_entries]
        elif kmn in Name_matched_KMN_KEGG_dict.keys():
            kegg_entries = [Name_matched_KMN_KEGG_dict[kmn]]
            kegg_names = [kegg_entry_name_dict[kegg] for kegg in kegg_entries]           

        else:
            kegg_entries = []
            kegg_names = []
        kmn_disease = KMN_all_df.at[kmn, "Disease_english"]

        if kmn in list(KMN_ICD10_dict.keys()):
            kmn_icd10 = KMN_ICD10_dict[kmn]
        else:
            kmn_icd10 = "NO_DATA"
        KMN_extention_df.at[kmn, "Keymolnet_disease"] = kmn_disease
        KMN_extention_df.at[kmn, "ICD10"] = kmn_icd10
        KMN_extention_df.at[kmn, "KEGG_ENTRY"] = list_to_txt(kegg_entries)
        KMN_extention_df.at[kmn, "KEGG_Name"] = list_to_txt(kegg_names)
    
    print(KMN_extention_df)
    # KMN_extention_df.to_csv("KMN_extention_210723.tsv", sep="\t")

    ## print keymolnet==False and DrugBank == True
    print(kegg_df[kegg_df["keymolnet_id"].isin(["NO_DATA"]) & ~kegg_df["DrugBank"].isin(["NO_DATA"])])
    KMN_f_DB_t_df = kegg_df[kegg_df["keymolnet_id"].isin(["NO_DATA"]) & ~kegg_df["DrugBank"].isin(["NO_DATA"])]
    KMN_f_DB_t_df.to_csv("KEGG_DISEASEp_KMNfalse140_0628.txt", sep="\t")


if __name__ == "__main__":

    kegg_f = "KEGG_DISEASEplus_OMIM_HMDB_210621.tsv"
    matching_f = "../matching_output/matching_no_apostropheS.csv"
    KMN_ICD10_f = "../KeyMolNet_ICD10_tables/KMN_with_ICD10code_0723.csv"
    KMN_all_f = "../KeyMolNet_ICD10_tables/KMN_wICDcode_0723.txt"
    
    main(kegg_f, matching_f, KMN_ICD10_f, KMN_all_f)
