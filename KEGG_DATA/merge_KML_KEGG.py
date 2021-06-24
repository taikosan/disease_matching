import pandas as pd
import sys, csv

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

def read_KMNgml(fname, delimiter):
    id_list = []
    name_list = []
    count_list = []
    molecule_list = []
    KMN_dict = {}

    with open(fname, "r") as f:
        for row in csv.reader(f, delimiter=delimiter, lineterminator="\n"):
            id_list.append(row[0])
            name_list.append(row[1])
            count_list.append(row[2])
            molecule_list.append(row[3:])
            KMN_dict[row[0]] = (row[1], row[2], row[3:])
    return id_list, KMN_dict

def merge_KMN_and_KEGG_molecule_list(kegg_uniprot_list, kegg_chEBI_list, kmn_molecule_list):
    # 1. add prefix "CHEBI:" to kegg_chEBI_list
    kegg_chEBI_list = [f"CHEBI:{x}" for x in kegg_chEBI_list]
    print(len(kegg_uniprot_list))
    print(len(kegg_chEBI_list))
    print(len(kmn_molecule_list))
    
    # 2. marge molecule list
    KMN_KEGG_molecule_list = list(set(kegg_uniprot_list + kegg_chEBI_list + kmn_molecule_list))
    print(len(KMN_KEGG_molecule_list))

    return KMN_KEGG_molecule_list
    
def main(kegg_file, kmn_file):
    # read KMN_file 
    KMN_id_list, KMN_dict = read_KMNgml(kmn_file, delimiter="\t")

    # read kegg file
    kegg_df = pd.read_table(kegg_file, index_col=0)
    print(kegg_df)
    # get KMN == True KEGG entries
    kegg_KMN_df = kegg_df[~kegg_df["keymolnet_id"].isin(["NO_DATA"])]
    print(kegg_KMN_df)

    ### merge KMN and KEGG molecule info

    with open("./merged_output/KeyMolNet_Processed_wKEGG_210624.txt", "w") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")

        countdata_df = pd.DataFrame(index=list(kegg_KMN_df.index))
        for kegg_idx in list(kegg_KMN_df.index):
            print(kegg_idx)
            ## KMN_id
            kmn_idx = kegg_KMN_df.at[kegg_idx, "keymolnet_id"]
            print(kmn_idx)

            # get kegg uniprot list
            kegg_uniprot_list = txt_to_list(kegg_KMN_df.at[kegg_idx, "UniProtID"])
            # get kegg(hmdb) chEBI_list
            kegg_chEBI_list = txt_to_list(kegg_KMN_df.at[kegg_idx, "HMDB_ChEBI"])

            ## get KMN molecule_list
            kmn_molecule_list = KMN_dict[kmn_idx][2]

            ## merge kegg and kml molecules
            kegg_kmn_molecule_list = merge_KMN_and_KEGG_molecule_list(kegg_uniprot_list, kegg_chEBI_list, kmn_molecule_list)


            ## add "KMN_KEGG_molecules" columns
            kegg_KMN_df.at[kegg_idx, "KEGG_KMN_count"] = len(kegg_kmn_molecule_list)
            kegg_KMN_df.at[kegg_idx, "KEGG_KMN_molecules"] = list_to_txt(kegg_kmn_molecule_list)
            
            ## add count columns to count_df
            countdata_df.at[kegg_idx, "KMN_id"] = kmn_idx
            countdata_df.at[kegg_idx, "Name"] = KMN_dict[kmn_idx][0]
            countdata_df.at[kegg_idx, "KeyMolNet"] = KMN_dict[kmn_idx][1]
            countdata_df.at[kegg_idx, "KEGG_Uniprot"] = kegg_KMN_df.at[kegg_idx, "UniProt_count"]
            countdata_df.at[kegg_idx, "HMDB_ChEBI"] = kegg_KMN_df.at[kegg_idx, "HMDB_ChEBI_count"]
            countdata_df.at[kegg_idx, "KEGG_KeyMolNet"] = len(kegg_kmn_molecule_list)
            ## write new_gml file
            writer.writerow([KMN_dict[kmn_idx][0], len(kegg_kmn_molecule_list)] + kegg_kmn_molecule_list)

        ## save kegg_KMN df and countdata_df
        print(kegg_KMN_df)
        print(countdata_df)
        countdata_df.index.name="KEGG_ENTRY"
        kegg_KMN_df.to_csv("./merged_output/KEGG_KMN_marged_210624.tsv", sep="\t")
        countdata_df.to_csv("./merged_output/KEGG_KMN_countdata_210624.tsv", sep="\t")
    
if __name__ == "__main__":

    kegg_file = "KEGG_DISEASEplus_OMIM_HMDB_KMN_210623.tsv"
    kmn_file = "../matching_output/KeyMolNet_processed_wID.txt"
    
    main(kegg_file, kmn_file)

