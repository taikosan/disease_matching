import pandas as pd
import sys

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

def read_HMDB_data(hmdb_fname):
    HMDB_df = pd.read_csv(hmdb_fname)
    ## remove ChEBI = NaN
    HMDB_df = HMDB_df.fillna("NO_DATA")
    HMDB_chEBI_df = HMDB_df[~HMDB_df["chebi_id"].isin(["NO_DATA"])]

    ## set MIMid as index, change datatype as "str"
    HMDB_chEBI_df = HMDB_chEBI_df.astype({"omim_id":str})
    HMDB_chEBI_df = HMDB_chEBI_df.reset_index().set_index("omim_id")
    HMDB_MIMlist = list(HMDB_chEBI_df.index)
    print(HMDB_MIMlist[0:10])
    return HMDB_chEBI_df, HMDB_MIMlist
    
def read_OMIMdb_data(OMIMdb_fname):
    OMIM_df = pd.read_excel("Chugai-QSP_OMIM_data_7118diseases.xlsx", index_col=0)
    ## remove "Uniprot-id(s)" =NaN
    OMIM_df = OMIM_df.fillna("NO_DATA")
    OMIM_uniprot_df = OMIM_df[~OMIM_df["Uniprot-id(s)"].isin(["NO_DATA"])]
    ## set Disease-MIMid as index, change data type as str
    OMIM_uniprot_df = OMIM_uniprot_df.astype({"Disease-MIMid":str})
    print(OMIM_uniprot_df["Disease-MIMid"].dtype)
    OMIM_uniprot_df = OMIM_uniprot_df.reset_index().set_index("Disease-MIMid")
    OMIM_db_MIMlist = list(OMIM_uniprot_df.index)
    print(OMIM_db_MIMlist[0:10])
    return OMIM_uniprot_df, OMIM_db_MIMlist
    
def read_KEGGdb_data(KEGGdb_fname):
    KEGG_df = pd.read_table(KEGGdb_fname, index_col=0)
    KEGG_omim_df = KEGG_df[~KEGG_df["OMIM"].isin(["NO_DATA"])]
    print(KEGG_omim_df)
    
    ### test df
    title_list = list(KEGG_omim_df.index)[0:20]
    test_df = KEGG_omim_df.loc[title_list, :]
    print(test_df)
    # return test_df
    return KEGG_df

def main(HMDB_fname, KEGdb_fname):
    HMDB_chEBI_df, HMDB_MIMlist = read_HMDB_data(HMDB_fname)
    KEGG_df = read_KEGGdb_data(KEGGdb_fname)    
    
    for idx in list(KEGG_df.index):
        print(idx, KEGG_df.at[idx, "Name"])
        HMDB_chEBI_list = []

        ## check "OMIM" == NO_DATA or not
        MIM_list = txt_to_list(KEGG_df.at[idx, "OMIM"])
        if MIM_list == "NO_DATA":
            continue
        else:
            print(MIM_list)
            ## create OMIM_list
            for mim in MIM_list:
                ### check if MIM_id is in HMDB_MIMlist or not
                if mim in HMDB_MIMlist:
                    ## add OMIM_list to mim_genes
                    tmp_HMDB_chEBI_list = txt_to_list(HMDB_chEBI_df.at[mim, "chebi_id"])
                    HMDB_chEBI_list.extend(tmp_HMDB_chEBI_list)
                else:
                    print("NOT in HMDB_chEBI list")
                    continue
            print("HMDB_mims:")
            print("HMDB_chEBI=", HMDB_chEBI_list)
            
        ## update KEGG_df
        if len(HMDB_chEBI_list) == 0:
            KEGG_df.at[idx, "HMDB_ChEBI"] = "NO_DATA"
            KEGG_df.at[idx, "HMDB_ChEBI_count"] = 0 
        else:
            KEGG_df.at[idx, "HMDB_ChEBI"] = list_to_txt(HMDB_chEBI_list)
            KEGG_df.at[idx, "HMDB_ChEBI_count"] = len(HMDB_chEBI_list)            

    ## save KEGG_df
    KEGG_df.to_csv("KEGG_DISEASEplus_OMIM_HMDB_210621.tsv", sep="\t")

if __name__ == "__main__":

    HMDB_fname = "hmdb_omim_acc.csv"
    KEGGdb_fname = "KEGG_DISEASEplus_OMIM_210617.tsv"

    main(HMDB_fname, KEGGdb_fname)
