import pandas as pd
import csv, sys

def read_gml(fname, delimiter):
    name_list = []
    count_list = []
    molecule_list = []
    with open(fname, "r") as f:
        for row in csv.reader(f, delimiter=delimiter, lineterminator="\n"):
            name_list.append(row[0])
            count_list.append(row[1])
            molecule_list.append(row[2:])
    return name_list, count_list, molecule_list

def match_ID_and_ProcessedName(processed_fname, ID_fname, new_gml_fname):
    ## get processed disease name list
    processed_disease_name_list, count_list, molecule_list = read_gml(processed_fname, delimiter="\t") 
    print(processed_disease_name_list[0:20])

    ## get ID and disease name(original name)
    id_df = pd.read_csv(ID_fname, index_col=0)
    print(id_df)

    ## replace white-space, "'" and "," to "_" on original name, and create replaced name - id dict
    name_id_dict = {}
    for idx in list(id_df.index):
        org_name = id_df.at[idx, "Disease_english"]
        replaced_name = org_name.translate(str.maketrans({" ":"_", "'":"_", ",":"_"}))
        # replaced_name = org_name.replace(" ", "_")
        # replaced_name = replaced_name.replace("'", "_")
        name_id_dict[replaced_name] = idx

    ## check processed disease names are in name_id_dict, and write new gml file with ID
    with open(new_gml_fname, "w") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for name, count, molecules in zip(processed_disease_name_list, count_list, molecule_list):
        # for name, count, molecules in zip(processed_disease_name_list[0:20], count_list[0:20], molecule_list[0:20]):
            if name == "Disease_english":
                print(name)
                print(count)
                print(molecules)
                continue
            KMN_id = name_id_dict[name]
            writer.writerow([f"KMN{KMN_id}", name, count] + molecules)
    
if __name__ == "__main__":
    
    processed_fname = "./disease_raw/KeyMolNet_processed.txt"
    ID_fname = "./disease_raw/Disease ID_KeyMolnet.csv"
    new_gml_fname = "./matching_output/KeyMolNet_processed_wID.txt"
    
    match_ID_and_ProcessedName(processed_fname, ID_fname, new_gml_fname)
