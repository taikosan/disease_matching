"""
Created on 2021/06/07

@author: akirakaren
"""
import pandas as pd
import os

from find_matching_noapostropheS import find_matching_noapostropheS
from find_matching_all_lowercase import find_matching_all_lowercase

### Import KEGG database and KeyMolnet database
KEGG = pd.read_csv("disease_raw/kegg_DISEASE_210514.csv")
KeyMolnet = pd.read_csv("disease_raw/Disease ID_keymolnet.csv")

### Make output directory
dirname = "./matching_output"
os.makedirs(dirname, exist_ok=True)

### Find matching diseases
matching_all_lowercase = find_matching_all_lowercase(KEGG, KeyMolnet,
                                       dirname, "matching_all_lowercase.csv")
matching_noapostorpheS = find_matching_noapostropheS(KEGG, KeyMolnet, dirname)

### Number of matching diseases based on change
with open(os.path.join(dirname, "metadata.txt"), "w+") as f:
    f.truncate(0)
    f.write("KEGG: 'kegg_DISEASE_210514.csv', " + str(len(KEGG)) + " diseases"
            "\nKeyMolnet: 'Disease ID_keymolnet.csv', " + str(len(KeyMolnet)) +
            " diseases  ")

    f.write("\n\n# of matching diseases: "
            "\nall_lowercase = " + str(len(matching_all_lowercase)) +
            "\nall_lowercase + no 's = " + str(len(matching_noapostorpheS)))

if __name__ == '__main__':
    pass
