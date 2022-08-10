

import pandas as pd
from functools import reduce
import json

# From config first, then from rule
output_folder = snakemake.config["out"] + snakemake.config["path_merged_data"]
sample_id = snakemake.config["sample_id"]

dataset_file_list = snakemake.input["dataset_list"]


# TODO metadata-pcas needs be clear as this can be done on other datasets
# TODO add name of dataset to folder and file names to stay organized

dataset_json = {
  "datasets": {
    "dataset": [
      {
        "ds_name": "DS_ASVS",
        "dim_method": "SCNIC",
        "file_path": "random-forest/TEST/HDL/01-DIM-SCNIC-biom/SCNIC_modules_for_workflow.txt",
        "param_dict": ""
      },
      {
        "ds_name": "DS_metadata",
        "dim_method": "PCA",
        "file_path": "data/hdl-test/real-data-no-asvs.txt",
        "param_dict": ""
      }
    ]
  }
}

dataset_json = json.dumps(dataset_json)
dataset_json = json.loads(dataset_json)

df_list = []
df_dict = {}

unique_cols = []
for dataset in dataset_json:
    for k, v in dataset_json[dataset].items():
        for i in v:
            df = pd.read_csv(i["file_path"], sep="\t",
                             index_col=sample_id)
            df_dict[i["ds_name"]] = df
            unique_cols += list(df.columns)
            df_list.append(df)

# TODO

if len(set(unique_cols)) == len(unique_cols):
    pass
else:
   print("Warning: Dataset columns are not unique and may create problems")

# TODO VERIFY THIS MERGE *** review merge decisions ***
df = reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True,
                                  how="inner", validate="one_to_one",
                                  sort=False), df_list)
#
# df = df[df.columns.drop(list(df.filter(regex='_reference')))]
print(df)

def _subset_simple(df_input, drop_rows, subset_rows, drop_cols, constrain_cols):

    for k, v in drop_rows.items():
        try:
            df_input = df_input.loc[(df_input[k] != v)]
        except:
            pass

    for k, v in subset_rows.items():
        try:
            df_input = df_input.loc[(df[k] == v)]
        except:
            pass

    if constrain_cols and drop_cols:
        print("Cannot have features (arguments) in 'drop_cols' and "
              "'constrain_cols' parameters")
        return

    if constrain_cols:
        df_input = df_input[constrain_cols]
    else:
        pass

    df_input = df_input.drop(drop_cols, axis=1)

    return df_input


drop_rows = {"SexualClassification": "Women"}
constrain_rows = {"Diet": "Agrarian"}
drop_cols = ["Inflammation", "TotalCholesterol", "PCA", "HOMAIR"]
constrain_cols = []

df = _subset_simple(df, drop_rows, constrain_rows, drop_cols,
                    constrain_cols)

df.to_csv(output_folder + "final-merged-dfs.txt", sep="\t")

