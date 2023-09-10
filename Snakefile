import json
import os
import sys

from snakemake.io import load_configfile
from snakemake import parse_config


# TODO improve these path numbers
config["path_dim_pca"] = "01-DIM-PCA/"
config["path_dim_scnic"] = "01-DIM-SCNIC/"

config["path_original"] = "04-SELECTED-FEATURES-original/"
config["path_first"] = "04-SELECTED-FEATURES-first/"
config["path_previous"] = "04-SELECTED-FEATURES-previous/"
config["path_pairwise"] = "04-SELECTED-FEATURES-pairwise/"

config["path_post_hoc"] = "05-POST-HOC-TESTS/"

config["version"] = "2023.08.0"

out_dir = config["out"]

# get the argument after the --configfile param
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

path_post_hoc = out_dir + config["path_post_hoc"]


def dim_reduction():
    dataset_json = json.loads(config["dataset_json"])
    file_path_list = []
    ds_param_dict_list = []
    pca_ds_list = []
    scnic_ds_list = []

    for ds_name in dataset_json["datasets"]:
        dim_method = dataset_json["datasets"][ds_name]["dim_method"]
        ds_param_dict_list.append(dataset_json["datasets"][ds_name]
        ["ds_param_dict"]["df_mod"])
        if dim_method == "pca":
            f_out = path_dim_pca + ds_name + "/PCA-dim-reduction-" \
                                             "for-workflow.txt"
            pca_ds_list.append(ds_name)

        elif dim_method == "scnic":
            f_out = path_dim_scnic + ds_name + \
                    "/SCNIC-dim-reduction-for-workflow.txt"
            scnic_ds_list.append(ds_name)

        else:
            f_out = dataset_json["datasets"][ds_name]["file_path"]

        file_path_list.append(f_out)

    return file_path_list, ds_param_dict_list, \
           pca_ds_list, \
           scnic_ds_list, dataset_json

config["file_path_list"], config["ds_param_dict_list"], \
config["pca_ds_list"], \
config["scnic_ds_list"], dataset_json = dim_reduction()
config["dataset_json"] = dataset_json


rule all:
    input:
        original = path_original + "original.pdf",
        first = path_first + "first.pdf",
        previous = path_previous + "previous.pdf",
        pairwise = path_pairwise + "pairwise.pdf"


def get_pca_in_file(wildcards):
    in_file = dataset_json["datasets"][wildcards.pca_ds_name]["file_path"]
    return in_file


def get_pca_groups_list(wildcards):
    pca_groups_list = dataset_json["datasets"][wildcards.pca_ds_name]["dim_param_dict"]["pca_groups_list"]
    return pca_groups_list


def get_scnic_in_file(wildcards):
    in_file = dataset_json["datasets"][wildcards.scnic_ds_name]["file_path"]
    return in_file


def get_scnic_method(wildcards):
    method = dataset_json["datasets"][wildcards.scnic_ds_name]["dim_param_dict"]["method"]
    return method


def get_scnic_min_r(wildcards):
    min_r = dataset_json["datasets"][wildcards.scnic_ds_name]["dim_param_dict"]["min_r"]
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
        pca_groups_list = get_pca_groups_list,
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
rule integrate_datasets:
    input:
        file_path_list = config["file_path_list"],
    output:
        out_file = path_original + "original.txt"
    params:
        ds_param_dict_list = config["ds_param_dict_list"],
        build_datatable = config["build_datatable"],
    conda: "conda_envs/r_env.yaml",
    script:
        "scripts/integrate_datasets.R"


# Create "delta datasets"; differences between samples from different
# timepoints, within individual
rule make_delta_datasets:
    input:
        in_file = path_original + "original.txt"
    output:
        out_file = out_dir + "04-SELECTED-FEATURES-{reference}/"
                   "{reference}.txt"
    params:
        reference_time = "{reference}",
        absolute_values = "no",
        build_datatable = config["build_datatable"], # TODO have this be optional default true
        distance_matrices = config["distance_matrices"],
    conda: "conda_envs/r_env.yaml",
    script:
        "scripts/create_deltas.R"


# Python in a complicated environment
# MERF has extremely specific Python requirements that could pose problems
# when combined with other software? Unsure.
# TODO add 'previous' or 'first' and other options to log file/report
# TODO prevent them from running pairwise from 'previous' file?
# TODO fix underscores in names.
#
# Also, this whole thing is done a little strange with the variables in the
# name of the output. Not sure how I want to handle this, but it was because
# I wanted my file names to reflect the arguments


# TODO check constrain, and wildcards (not a good use)
rule random_forest:
    input:
        in_file = out_dir + "04-SELECTED-FEATURES-{reference}/{reference}.txt"
    output:
        out_file = out_dir + "04-SELECTED-FEATURES-{reference}/{reference}.pdf"
    params:
        dataset = "{reference}",
        random_effect = config["random_effect"],
        sample_id = config["sample_id"],
        response_var = config["response_var"],
        include_time = config["include_time"],
        max_features = config["max_features"],
        n_estimators = config["n_estimators"],
        iterations = config["iterations"],  # 20 suggested, 10 testing (MERF)
        borutashap_threshold = config["borutashap_threshold"],
        borutashap_p = config["borutashap_p"],
        borutashap_trials = config["borutashap_trials"]
    conda: "conda_envs/merf.yaml"
    script:
        "scripts/random_forest.py"

# TODO see if I can create these files for deltas on original timepoint check
# feature summaries
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
        "scripts/feature-heatmaps.R"


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
        "scripts/url_interpretation.py"

# TODO program if statements here and make report reflect missing data
# TODO consider if statements with all rule
# TODO idea make four new render reports
# TODO consider making some of these parameters
# rule run_post_hoc_stats:
#     input:
#         deltas_first = path_first + "first.txt",
#         deltas_previous = path_previous + "previous.txt",
#         deltas_pairwise = path_pairwise + "pairwise.txt",
#         df_original = path_original + "original.txt",
#         fixed_effects_first = path_first + "first-boruta-important.txt",
#         fixed_effects_previous = path_previous + "previous-boruta-important.txt",
#         fixed_effects_pairwise = path_pairwise + "pairwise-boruta-important.txt",
#         fixed_effects_original = path_original + "original-boruta-important.txt"
#     output:
#         out_file = path_post_hoc + "post-hoc-analysis.html"
#     params:
#         response_var = config["response_var"],
#     conda: "conda_envs/r_env.yaml",
#     script:
#         "scripts/post_hoc_tests.R"

# Only require log files because then report will always render
# TODO add error reports to log file so people know entire 
# algorithm did not work
# TODO remove feature occurrances possibly
rule render_report:
    input:
        original_log = path_original + "original-log.txt",
        first_log = path_first + "first-log.txt",
        previous_log= path_previous + "previous-log.txt",
        pairwise_log = path_pairwise + "pairwise-log.txt",
        feature_occurances = out_dir + "important-feature-occurrences.svg",
        urls = out_dir + "urls.txt",
        # post_hoc = path_post_hoc + "post-hoc-analysis.html",
        # post_hoc_viz = path_post_hoc + "post-hoc-combined.pdf"
    output:
        md_doc = out_dir + "EXPLANA-report.html",
    conda: "conda_envs/r_env.yaml",
    script:
        "scripts/report.Rmd"