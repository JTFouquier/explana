import json
import os
import sys
import shutil

from snakemake.io import load_configfile
from snakemake import parse_config

if not os.path.exists(config["out"]):
    os.makedirs(config["out"])

index = sys.argv.index("--configfile")
config_file = sys.argv[index + 1]
out_dir = config["out"]

# make a copy of the config file for repository
new_filename = "analysis_config.yaml" 
shutil.copy(config_file, os.path.join(out_dir, new_filename))

config["path_dim_pca"] = "DIM-PCA/"
config["path_dim_scnic"] = "DIM-SCNIC/"

# can simplifity this code
config["path_original"] = "SELECTED-FEATURES-original/"
config["path_first"] = "SELECTED-FEATURES-first/"
config["path_previous"] = "SELECTED-FEATURES-previous/"
config["path_pairwise"] = "SELECTED-FEATURES-pairwise/"

config["version"] = "2023.11.0"


rule all:
    input:
        original = f"{out_dir}SELECTED-FEATURES-original/original.pdf",
        first = f"{out_dir}SELECTED-FEATURES-first/first.pdf",
        previous = f"{out_dir}SELECTED-FEATURES-previous/previous.pdf",
        pairwise = f"{out_dir}SELECTED-FEATURES-pairwise/pairwise.pdf",
        md_doc = f"{out_dir}EXPLANA-report.html"


index = sys.argv.index("--configfile")
config_file = sys.argv[index + 1]
cl_configfile = dict(load_configfile(config_file))

config.update(cl_configfile)

# paths with the final workflow results path included
path_original = out_dir + config["path_original"]
path_first = out_dir + config["path_first"]
path_previous = out_dir + config["path_previous"]
path_pairwise = out_dir + config["path_pairwise"]

path_dim_pca = out_dir + config["path_dim_pca"]
path_dim_scnic = out_dir + config["path_dim_scnic"]


def dim_reduction():
    input_datasets = config["input_datasets"]
    file_path_list = []
    ds_param_dict_list = []
    pca_ds_list = []
    scnic_ds_list = []

    for ds_name in input_datasets:
        dim_method = input_datasets[ds_name]["dim_method"]
        ds_param_dict_list.append(input_datasets[ds_name]["df_mod"])
        if dim_method == "pca" or dim_method == "PCA":
            f_out = path_dim_pca + ds_name + "/PCA-dim-reduction-" \
                                             "for-workflow.txt"
            pca_ds_list.append(ds_name)

        elif dim_method == "scnic" or dim_method == "SCNIC":
            f_out = path_dim_scnic + ds_name + \
                    "/SCNIC-dim-reduction-for-workflow.txt"
            scnic_ds_list.append(ds_name)

        elif dim_method == "None":
            f_out = input_datasets[ds_name]["file_path"]

        file_path_list.append(f_out)

    return file_path_list, ds_param_dict_list, \
           pca_ds_list, \
           scnic_ds_list, input_datasets

config["file_path_list"], config["ds_param_dict_list"], \
config["pca_ds_list"], \
config["scnic_ds_list"], input_datasets = dim_reduction()
config["input_datasets"] = input_datasets


def get_pca_in_file(wildcards):
    in_file = input_datasets[wildcards.pca_ds_name]["file_path"]
    return in_file


def get_pca_list(wildcards):
    pca_list = input_datasets[wildcards.pca_ds_name]["dim_param_dict"]["pca_list"]
    return pca_list

def get_scnic_in_file(wildcards):
    in_file = input_datasets[wildcards.scnic_ds_name]["file_path"]
    return in_file

def get_scnic_method(wildcards):
    method = input_datasets[wildcards.scnic_ds_name]["dim_param_dict"]["method"]
    return method


def get_scnic_min_r(wildcards):
    min_r = input_datasets[wildcards.scnic_ds_name]["dim_param_dict"]["min_r"]
    return min_r


# Perform dimensionality reduction using principal coordinates analysis (PCA)
# Can perform PCA on as many variables as you would like using config file
rule dim_reduction_pca:
    input:
        in_file = get_pca_in_file,
    output:
        out_file = path_dim_pca + "{pca_ds_name}" + "/PCA-dim-reduction-"
                                                    "for-workflow.txt"
    params:
        dataset_name = expand("{pca_ds_name}",
            pca_ds_name=config["pca_ds_list"]),
        pca_list = get_pca_list,
    conda: "conda_envs/r_env.yaml",
    script:
        "scripts/dim_reduction_pca.R"


# TODO remove in_file_metadata
rule dim_reduction_scnic:
    input:
        in_file = get_scnic_in_file,
    output:
        out_file =  path_dim_scnic + "{scnic_ds_name}" + "/SCNIC-dim-reduction"
                                                         "-for-workflow.txt"
    params:
        dataset_name = expand("{scnic_ds_name}",
            scnic_ds_name=config["scnic_ds_list"]),
        method = get_scnic_method,
        min_r = get_scnic_min_r
    conda: "conda_envs/scnic.yaml",
    script:
        "scripts/dim_reduction_scnic.py"


# After doing dimensinoality reductions and filters on individual datasets
# TODO handle biom tables; for now, must be 
rule integrate_datasets:
    input:
        file_path_list = config["file_path_list"],
    output:
        out_file = path_original + "original.txt"
    params:
        ds_param_dict_list = config["ds_param_dict_list"],
    conda: "conda_envs/r_env.yaml",
    script:
        "scripts/integrate_datasets.R"


# Create "delta datasets"; differences between samples from different
# timepoints, within individual
rule make_delta_datasets:
    input:
        in_file = path_original + "original.txt"
    output:
        out_file = out_dir + "SELECTED-FEATURES-{reference}/"
                   "{reference}.txt"
    params:
        reference_time = "{reference}",
        include_reference_values = config["include_reference_values"],
        absolute_values = config["absolute_values"],
        distance_matrices = config["distance_matrices"],
    conda: "conda_envs/r_env.yaml",
    script:
        "scripts/create_deltas.R"


rule random_forest:
    input:
        in_file = out_dir + "SELECTED-FEATURES-{reference}/{reference}.txt"
    output:
        out_file = out_dir + "SELECTED-FEATURES-{reference}/{reference}.pdf",
        out_boruta = out_dir + \
        "SELECTED-FEATURES-{reference}/{reference}-boruta-important.txt",
    params:
        dataset = "{reference}",
        random_effect = config["random_effect"],
        sample_id = config["sample_id"],
        response_var = config["response_var"],
        include_time = config["include_time"],
        max_features = config["max_features"],
        n_estimators = config["n_estimators"],
        iterations = config["iterations"],
        borutashap_threshold = config["borutashap_threshold"],
        borutashap_p = config["borutashap_p"],
        borutashap_trials = config["borutashap_trials"]
    conda: "conda_envs/merf.yaml"
    script:
        "scripts/random_forest.py"


rule final_steps_r:
    input:
        selected_features_original = path_original + "original-boruta-important.txt",
        selected_features_first = path_first + "first-boruta-important.txt",
        selected_features_previous = path_previous + "previous-boruta-important.txt",
        selected_features_pairwise = path_pairwise + "pairwise-boruta-important.txt"
    output:
        out_file = out_dir + "all-important-features.txt",
        out_file_image = out_dir + "important-feature-occurrences.svg",
    conda: "conda_envs/r_env.yaml",
    script:
        "scripts/final_steps.R"


rule final_steps_python:
    input: 
        feature_file = out_dir + "all-important-features.txt",
        original_features = path_original + "original-boruta-important.txt",
        first_features = path_first + "first-boruta-important.txt",
        previous_features = path_previous + "previous-boruta-important.txt",
        pairwise_features = path_pairwise + "pairwise-boruta-important.txt"
    output:
        out_file = out_dir + "urls.txt",
    conda: "conda_envs/merf.yaml",
    script:
        "scripts/final_steps.py"


rule render_report:
    input:
        original_pdf = path_original + "original.pdf",
        first_pdf = path_first + "first.pdf",
        previous_pdf = path_previous + "previous.pdf",
        pairwise_pdf = path_pairwise + "pairwise.pdf",
        feature_occurances = out_dir + "important-feature-occurrences.svg",
        urls = out_dir + "urls.txt",
    output:
        md_doc = out_dir + "EXPLANA-report.html",
    conda: "conda_envs/r_env.yaml",
    script:
        "scripts/report.Rmd"