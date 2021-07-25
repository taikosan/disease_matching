import pandas as pd
import csv, sys

def read_gml(fname, delimiter):
    id_list = []
    name_list = []
    count_list = []
    molecule_list = []
    with open(fname, "r") as f:
        for row in csv.reader(f, delimiter=delimiter, lineterminator="\n"):
            id_list.append(row[0])
            name_list.append(row[1])
            count_list.append(row[2])
            molecule_list.append(row[3:])
    return id_list, name_list, count_list, molecule_list

def create_ICDcode_dict(df):
    print(df)
    name_id_dict = {}
    for idx in list(df.index):
        ## 1. replace all letters to lowercase
        orgname = df.at[idx, "Name1"]
        orgname = orgname.lower()
        ## 2. replace white-space, "'" and "," to "_" on original name, and create replaced name - id dict
        replaced_name = orgname.translate(str.maketrans({" ":"_", "'":"_", ",":"_"}))
        name_id_dict[replaced_name] = df.at[idx, "icd10_code"]

    return name_id_dict

def main(icd10_f, KMN_f):

    ## read icd10 code and create icdcode-name dict
    icd10_df = pd.read_csv(icd10_f, index_col=0)
    print(icd10_df)
    icd_name_dict= create_ICDcode_dict(icd10_df)

    print(list(icd_name_dict.keys())[10:20])

    ### check names
    with open("test.csv", "w") as out:
        writer = csv.writer(out)
        writer.writerow(["Name","ICD10"])
        for name in list(icd_name_dict.keys()):
            writer.writerow([name, icd_name_dict[name]])

    ## read KMN file
    KMN_IDs, KMN_names, KMN_counts, KMN_items = read_gml(KMN_f, "\t")
    print(KMN_names[0:10])
    
    ## convert KMN names to lowercase
    KMN_lower_names = [x.lower() for x in KMN_names]

    ## generate KMN-icd10code 
    KMN_icd10codes = []
    for KMNname in KMN_lower_names:
        print(KMNname)
        if KMNname in list(icd_name_dict.keys()):
            KMN_icd10codes.append(icd_name_dict[KMNname])
        else:
            KMN_icd10codes.append("NO_DATA")


    with open("./matching_output/KeyMolNet_processed_icd10.txt" ,"w") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for _id, name, icd10, count, molecules in zip(KMN_IDs, KMN_names, KMN_icd10codes, KMN_counts, KMN_items):
            writer.writerow([_id, name, icd10, count] + molecules)

    with open("./matching_output/KeyMolNet_icd10_noItems.txt" ,"w") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["KMN_ID", "Name", "ICD10", "Count"])
        for _id, name, icd10, count in zip(KMN_IDs, KMN_names, KMN_icd10codes, KMN_counts):
            writer.writerow([_id, name, icd10, count])    
    
if __name__ == "__main__":

    icd10_f = "./ICD10_code/ICD10_code_Name_2021.csv" ## dataframe
    KMN_f = "./matching_output/KeyMolNet_processed_wID.txt" ## gml file

    main(icd10_f, KMN_f)
