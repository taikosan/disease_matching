"""
Created on 2021/06/07

@author: akirakaren
"""
import pandas as pd

from find_matching_all_lowercase import find_matching_all_lowercase


def find_matching_noapostropheS(kegg, keymolnet, dirname):
    """Remove apostrophe s from KeyMolnet's disease names.
    Run find_matching_all_lowercase()"""
    keymolnet["Disease_english"] = keymolnet["Disease_english"].str.replace(
        "'s", "")
    return find_matching_all_lowercase(kegg, keymolnet, dirname, "matching_no_apostropheS.csv")


if __name__ == '__main__':
    pass
