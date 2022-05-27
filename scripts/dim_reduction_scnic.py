

import os
from subprocess import Popen, PIPE
import time

import pandas as pd
from biom import load_table

# TODO fix this: prefix + process + dataset as folder
scnic_out_folder = "random-forest/TEST/HDL/01-DIM-SCNIC--biom/"
# TODO rename
biom_name = snakemake.input["in_file"]
in_file_metadata = snakemake.input["in_file_metadata"]


def makeSafeDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def within(out_folder, biom_file_name):
    process = Popen(["SCNIC_analysis.py", "within", "-i", biom_file_name,
                     "-o", out_folder + "within_output/", "-m", "spearman"],
                    stdout=PIPE, stderr=PIPE)
    process.communicate()


def modules(out_folder, biom_file_name):
    process = Popen(["SCNIC_analysis.py", "modules", "-i", out_folder +
                     "within_output/correls.txt", "-o", out_folder +
                     "modules_output/", "--min_r", "0.45", "--table",
                     biom_file_name], stdout=PIPE, stderr=PIPE)
    process.communicate()


def main(scnic_folder, biom_file_name):
    makeSafeDir(scnic_folder)

    # run SCNIC
    within(scnic_folder, biom_file_name)
    modules(scnic_folder, biom_file_name)

    # load biom table, convert samples to rows; OTUs/ASVs as columns (for
    # adding to large df convert to df and add SCNIC modules to dataframe
    biom_table = load_table(scnic_folder + "modules_output/collapsed.biom")
    biom_table = biom_table.transpose()
    df_scnic = biom_table.to_dataframe().rename_axis("StudyID.Timepoint")
    df_scnic.to_csv(scnic_folder + "SCNIC_modules.txt", sep="\t")

    df_metadata = pd.read_csv(in_file_metadata, sep="\t",
                              index_col="SampleID")

    # TODO update sample naming system for workflow
    # TODO check params here; suffixes, duplicate columns, etc
    df = pd.merge(df_metadata, df_scnic, left_index=True, right_index=True,
                  how="inner", validate="one_to_one", sort=False)
    df = df.rename_axis("StudyID.Timepoint")
    df.to_csv(scnic_folder + "metadata_with_SCNIC_modules.txt", sep="\t")

start_time = time.time()
main(scnic_folder=scnic_out_folder, biom_file_name=biom_name)
print("--- %s seconds ---" % (time.time() - start_time))
print((time.time() - start_time)/60)












