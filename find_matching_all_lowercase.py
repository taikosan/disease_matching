"""
Created on 2021/06/07

@author: akirakaren
"""
import pandas as pd
import os


def disease_to_id(df, disease_name, disease_col_str, id_col_str):
    """Find disease ID from disease name."""
    idx = df.index[df[disease_col_str] == disease_name].tolist()[0]
    return df.iloc[idx, df.columns.get_loc(id_col_str)]


def find_matching_all_lowercase(kegg, keymolnet, dirname, outfile_str):
    """Find diseases that match between KEGG database and KeyMolnet database
    when disease name are all lowercase."""

    ### Convert disease name to lower case
    kegg['Name'] = kegg['Name'].str.lower()
    keymolnet["Disease_english"] = keymolnet["Disease_english"].str.lower()

    ### Find in common diseases
    incommon = list(set(kegg['Name']) & set(keymolnet["Disease_english"]))

    ### Find KEGG_ID and KeyMolnet_ID for each disease
    keggidlst, molnetidlst = [], []
    for name in incommon:
        keggid = disease_to_id(kegg, name, 'Name', 'ENTRY')
        keggidlst.append(keggid)

        molnetid = disease_to_id(keymolnet, name, 'Disease_english', 'No')
        molnetidlst.append(molnetid)

    ### Create Dataframe
    matchingdict = {'Disease_name': incommon, 'KEGG_ID': keggidlst,
                    'KeyMolnet_ID': molnetidlst}
    matching = pd.DataFrame(matchingdict)
    sorted_matching = matching.sort_values(by=['KEGG_ID'])

    ### Output to CSV
    # dirname = "./matching_output"
    # os.makedirs(dirname, exist_ok=True)
    filename = os.path.join(dirname, outfile_str)
    sorted_matching.to_csv(filename, index=False)
    return sorted_matching


if __name__ == '__main__':
    pass
    # ### Import KEGG database and KeyMolnet database
    # KEGG = pd.read_csv("disease_raw/KEGG_DISEASE_210514.csv")
    # KeyMolnet = pd.read_csv("disease_raw/Disease ID_KeyMolnet.csv")
    #
    # print(find_matching_all_lowercase(KEGG, KeyMolnet,
    #                                   "matching_all_lowercase.csv"))
