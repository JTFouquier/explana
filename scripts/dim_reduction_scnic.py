

import os
from subprocess import Popen, PIPE
import time

import pandas as pd
from biom import load_table

# TODO fix this: prefix + process + dataset as folder

scnic_out_folder = snakemake.config["out"] + \
                   snakemake.config["path_dim_scnic"] + \
                   snakemake.config["ds_name"] + "/"
biom_name = snakemake.input["in_file"]
sample_id = snakemake.config["sample_id"]


def makeSafeDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def within(out_folder, biom_file_name):
    process = Popen(["SCNIC_analysis.py", "within", "-i", biom_file_name,
                     "-o", out_folder + "within_output/", "-m", "spearman"],
                    stdout=PIPE, stderr=PIPE)
    process.communicate()


# TODO add params
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
    df_scnic = biom_table.to_dataframe().rename_axis(sample_id)
    df_scnic.to_csv(scnic_folder + "SCNIC_modules_for_workflow.txt", sep="\t")

start_time = time.time()
main(scnic_folder=scnic_out_folder, biom_file_name=biom_name)
print("--- %s seconds ---" % (time.time() - start_time))
print((time.time() - start_time)/60)












