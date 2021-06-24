import pandas as pd
import sys

def main(kegg_f, matching_f):
    ## read kegg file
    kegg_df = pd.read_table(kegg_f, index_col=0)
    ## add kegg df to "keymolnet_id" column
    kegg_df["keymolnet_id"] = "NO_DATA"
    print(kegg_df)

    ## read id matching file
    matching_df = pd.read_csv(matching_f, index_col="KEGG_ID")
    print(matching_df)
    
    ## check id matching 
    for idx in list(matching_df.index):
        keymolnet_id = matching_df.at[idx, "KeyMolnet_ID"]
        kegg_df.at[idx, "keymolnet_id"] = f"KMN{keymolnet_id}"
    print(kegg_df)
    
    ## print keymolnet == True
    print(kegg_df[~kegg_df["keymolnet_id"].isin(["NO_DATA"])])

    ## print keymolnet == True and DrugBank == True
    print(kegg_df[~kegg_df["keymolnet_id"].isin(["NO_DATA"]) & ~kegg_df["DrugBank"].isin(["NO_DATA"])])

    ## save kegg_df
    kegg_df.to_csv("KEGG_DISEASEplus_OMIM_HMDB_KMN_210623.tsv", sep="\t")
    


if __name__ == "__main__":

    kegg_f = "KEGG_DISEASEplus_OMIM_HMDB_210621.tsv"
    matching_f = "./disease_matching/matching_output/matching_no_apostropheS.csv"
    
    main(kegg_f, matching_f)
