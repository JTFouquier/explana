

import pandas as pd
from functools import reduce

# From config first, then from rule
output_folder = snakemake.config["out"] + snakemake.config["path_merged_data"]
sample_id = snakemake.config["sample_id"]

fp_list = snakemake.input["file_path_list"]


# TODO metadata-pcas needs be clear as this can be done on other datasets
# TODO add name of dataset to folder and file names to stay organized

df_list = []
df_dict = {}

unique_cols = []

for i in range(len(fp_list)):
    df = pd.read_csv(fp_list[i], sep="\t", index_col=sample_id)
    unique_cols += list(df.columns)
    df_list.append(df)

# TODO create fail here
if len(set(unique_cols)) == len(unique_cols):
    pass
else:
   print("Warning: Dataset columns are not unique and may create problems")

# TODO VERIFY THIS MERGE *** review merge decisions ***
# TODO improve memory use here, don't load all to df
df = reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True,
                                  how="inner", validate="one_to_one",
                                  sort=False), df_list)
#
# df = df[df.columns.drop(list(df.filter(regex='_reference')))]

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


drop_rows = eval(snakemake.config["drop_rows"])
constrain_rows = eval(snakemake.config["constrain_rows"])
drop_cols = eval(snakemake.config["drop_cols"])
constrain_cols = eval(snakemake.config["constrain_cols"])

df = _subset_simple(df, drop_rows, constrain_rows, drop_cols,
                    constrain_cols)

df.to_csv(output_folder + "final-merged-dfs.txt", sep="\t")

