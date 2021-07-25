import pandas as pd
from collections import defaultdict
import sys, pprint

def list_to_txt(_list):
    if len(_list) == 0:
        txt = "NO_DATA"
    else:
        txt = ";".join(_list)
    return txt

def txt_to_list(txt):
    if txt == "NO_DATA":
        _list = []
    else:
        _list = txt.split(";")
    ## remove NA
    _list = [x for x in _list if not x == "NA"]
    _list = [x for x in _list if not x == "nan"]
    return _list

def read_KMN_ICDlist(fname, kegg_kmn_fname):
    df = pd.read_csv(fname, index_col=0)
    print(df)

    ## get KEGG=True KMN id list
    kegg_df = pd.read_table(kegg_kmn_fname)
    KEGGtrue_id_list = kegg_df["keymolnet_id"]
    KEGGtrue_id_list = [int(x.lstrip("KMN")) for x in KEGGtrue_id_list]
    
    ## generate ICD==key, No,Name==value dict
    preKMN_dict = defaultdict(list)
    for idx in df.index:
        ## remove KEGG=True diseases
        if idx in KEGGtrue_id_list:
            continue
        else:
            preKMN_dict[df.at[idx, "ICD___"]].append((f"KMN{idx}", df.at[idx, "Disease_english"]))

    print(len(preKMN_dict.keys()))

    ## convert {key: [[id_list], [name_list]]}

    KMN_dict = {}
    for key in preKMN_dict.keys():
        id_list = []
        name_list = []
        for values in preKMN_dict[key]:
            id_list.append(values[0])
            name_list.append(values[1])
        KMN_dict[key] = (id_list, name_list)

    return KMN_dict

def list_item_filter(_list, _code):
    filt_items = []
    print(_code[0:2])
    for item in _list:
        if item.startswith(_code[0:2]):
            print(item)
            filt_items.append(item)
    print(filt_items)
    return filt_items

def split_icd_midclass(key):
    [key1, key2] = key.split("-")
    header = key1[0]
    f1 = float(key1[1:])
    f2 = float(key2[1:])
    return header, f1, f2
    
def main(kmn_icd_fname, kmn_kegg_merged_fname, kegg_140_diseases_fname):
    ## read KMN ICD list and  remove KEGG_true KMN IDs
    KMN_dict = read_KMN_ICDlist(kmn_icd_fname, kmn_kegg_merged_fname)
    KMN_dict_keys = list(KMN_dict.keys())
    pprint.pprint(KMN_dict_keys[0:10])
    
    ## read kegg_140 diseases list
    kegg_df = pd.read_table(kegg_140_diseases_fname, index_col=0)
    print(kegg_df["ICD-10"])
    
    tidy_list = []
    ###
    ## create tidy table of kegg_list
    for kegg_idx in kegg_df.index:
        icd_list = txt_to_list(kegg_df.at[kegg_idx, "ICD-10"])
        for icd in icd_list:
            tidy_list.append([kegg_idx, kegg_df.at[kegg_idx, "Name"], kegg_df.at[kegg_idx, "Alias_Name"], icd, icd[0], icd[1:]])
    # pprint.pprint(tidy_list)
    tidy_df = pd.DataFrame(tidy_list, columns=["ENTRY","KEGG_disease_name",  "KEGG_synosyms", "ICD-10", "header", "float"])
    tidy_df = tidy_df.astype({"float":float})
    print(tidy_df)

    tidy_df["ICD_mid_class"] = "NO_DATA"
    tidy_df["KMN_candidate_IDs"] = "NO_DATA"
    tidy_df["KMN_candidate_Names"] = "NO_DATA"

    # pprint.pprint(list(KMN_dict.values()))
    
    ## add tidy_df to KMN_dict_keys
    new_tidy_list = []
    for key in KMN_dict_keys:
        # key='C81-C96'
        header, f1, f2 = split_icd_midclass(key)
        header_df = tidy_df[tidy_df["header"].isin([header])]
        hit_index_list = list(header_df.query(f"{f1} <= float <= {f2}").index)
        print(tidy_df.loc[hit_index_list,:])
        for idx in hit_index_list:
            s = tidy_df.loc[idx]
            for kmn_id, kmn_name in zip(KMN_dict[key][0], KMN_dict[key][1]):
                s_cp = s.copy()
                s_cp["ICD_mid_class"] = key
                s_cp["KMN_candidate_IDs"] = kmn_id
                s_cp["KMN_candidate_Names"] = kmn_name
                new_tidy_list.append(s_cp)

    print(len(new_tidy_list))

    new_tidy_df = pd.concat(new_tidy_list, axis=1).T
    print(new_tidy_df)


    ## read KMN-ICD10 table
    KMN_ICD10df = pd.read_table("../matching_output/KeyMolNet_icd10_noItems.txt", index_col=0)
    print(KMN_ICD10df)
    
    new_tidy_df["KMN_ICD10"] = [ KMN_ICD10df.at[kmn_id, "ICD10"] for kmn_id in new_tidy_df["KMN_candidate_IDs"]]

    print(new_tidy_df)


    # new_tidy_df.to_csv("KMN_candidate_140_tidy_wICD10.txt", sep="\t")


    ### get unique KMN_candidate list
    KMN_candidate_table = new_tidy_df.loc[:, ["KMN_candidate_IDs", "KMN_candidate_Names", "ICD_mid_class", "KMN_ICD10"]]
    print(KMN_candidate_table)
    unique_candidate_table = KMN_candidate_table.drop_duplicates()
    print(unique_candidate_table)
    unique_candidate_table.columns = ["KMN_candidate_IDs", "KMN_candidate_Names", "ICDclass", "KMN_ICD10"]
    print(unique_candidate_table)
    
    # unique_candidate_table.to_csv("KMN_candidate_list.txt", "\t", index=False)
    
    sys.exit()
    
if __name__ == "__main__":

    kmn_icd_fname = "../disease_raw/Disease ID_KeyMolnet.csv"
    kmn_kegg_merged_fname = "../KEGG_DATA/merged_output/KEGG_KMN_marged_210624.tsv"

    kegg_140_diseases_fname = "KEGG_DISEASEp_KMNfalse140_0628.txt"
    
    
    main(kmn_icd_fname, kmn_kegg_merged_fname, kegg_140_diseases_fname)
