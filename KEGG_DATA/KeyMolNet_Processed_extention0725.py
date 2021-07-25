import pandas as pd
import sys, csv
import pprint

def get_KEGGdf_elements(kegg_df, kegg_entries_list, kegg_colname):
    target_element_list = []
    for kegg_entry in kegg_entries_list:
        target_list = txt_to_list(kegg_df.at[kegg_entry, kegg_colname])
        for elem in target_list:
            target_element_list.append(elem)
    ## remove duplicated elements
    target_element_list = sorted(list(set(target_element_list)))
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


def main(kegg_f, kmn_extention_f, kmn_processed_f):
    ## read kegg_f
    kegg_df = pd.read_table(kegg_f, index_col=0)
    ## read kmn_extention_f
    kmn_extention_df = pd.read_table(kmn_extention_f, index_col=0)
    print(kegg_df)
    print(kmn_extention_df)

    kegg_colname_list = ["MeSH", "OMIM", "Gene_symbol", "UniProtID", "Drug_name", "kegg_drugID", "DrugBank", "HMDB_ChEBI"]
    kmn_KeggElements_df = pd.DataFrame(data="NO_DATA", index=kmn_extention_df.index, columns=kegg_colname_list)
    ## add kmn disease name to kmn_KeggElements_df
    kmn_KeggElements_df["Keymolnet_disease"] = kmn_extention_df["Keymolnet_disease"]
    kmn_KeggElements_df = kmn_KeggElements_df[["Keymolnet_disease", "MeSH", "OMIM", "Gene_symbol", "UniProtID", "Drug_name", "kegg_drugID", "DrugBank", "HMDB_ChEBI"]]
    print(kmn_KeggElements_df)
    
    for kmn in list(kmn_extention_df.index):
        ## get KEGG entrie list
        kegg_entries_list = txt_to_list(kmn_extention_df.at[kmn, "KEGG_ENTRY"])
        for kegg_colname in kegg_colname_list:
            target_element = get_KEGGdf_elements(kegg_df, kegg_entries_list, kegg_colname)
            ## update kmn_KeggElements columns = kegg_colname
            kmn_KeggElements_df.at[kmn, kegg_colname] = target_element

    print(kmn_KeggElements_df)
    
    ## save kmn_KeggElements_df
    # kmn_KeggElements_df.to_csv("KMN_KeggElements_210725.tsv", sep="\t")


    ### create extended KeyMolNet_processed list
    with open(kmn_processed_f, "r") as f, open("KeyMolNet_processed_Extention_210725.txt", "w") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for row in csv.reader(f, delimiter="\t", lineterminator="\n"):
            print(row)
            kmn = row[0]
            elem_list = row[3:]
            kegg_UniProt_list = txt_to_list(kmn_KeggElements_df.at[kmn, "UniProtID"])
            kegg_ChEBI_list = txt_to_list(kmn_KeggElements_df.at[kmn, "HMDB_ChEBI"])
            ## add "ChEBI: to kegg_ChEBI elements"
            kegg_ChEBI_list = [f"CHEBI:{x}" for x in kegg_ChEBI_list]
            ## add KeggElements to KeyMolNet elements
            elem_list += kegg_UniProt_list
            elem_list += kegg_ChEBI_list
            ## remove duplicated elements
            elem_list = sorted(list(set(elem_list)))
            ## write new_row
            new_row = [kmn, row[1], len(elem_list)] + elem_list
            print(new_row)
            writer.writerow(new_row)


if __name__ == "__main__":
    kegg_f = "KEGG_DISEASEplus_OMIM_HMDB_KMN_210723.tsv"
    kmn_extention_f = "KMN_extention_210723.tsv"
    kmn_processed_f = "KeyMolNet_processed_wID.txt"

    main(kegg_f, kmn_extention_f, kmn_processed_f)
