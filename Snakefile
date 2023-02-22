import json
import os

# first configfile is the only thing you should change in snakefile;
# (note: data from both 'configfile' .yamls are added to the config)
configfile: "config/config.yaml"
# Do not modify paths or bad things will happen. :)
configfile: "config/config-paths.yaml"

# TODO add outputs concatenated with final output file
path_rf_original = config["out"] + config["path_rf_original"]
path_rf_first = config["out"] + config["path_rf_first"]
path_rf_previous = config["out"] + config["path_rf_previous"]
path_rf_pairwise = config["out"] + config["path_rf_pairwise"]

path_dim_pca = config["out"] + config["path_dim_pca"]
path_dim_scnic = config["out"] + config["path_dim_scnic"]

path_merged_data = config["out"] + config["path_merged_data"]

path_post_hoc = config["out"] + config["path_post_hoc"]

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
        original=path_rf_original + "MERF-original.pdf",
        first=path_rf_first + "mixed-RF-deltas-first.pdf",
        previous=path_rf_previous +
                 "mixed-RF-deltas-previous.pdf",
        pairwise=path_rf_pairwise +
                 "mixed-RF-deltas-pairwise.pdf"

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
        pca_groups_list = get_pca_groups_list
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
    conda: "conda_envs/scnic.yaml"
    script:
        "scripts/dim_reduction_scnic.py"

# After doing dimensinoality reductions and filters on individual datasets
rule integrate_datasets:
    input:
        file_path_list = config["file_path_list"],
    output:
        out_file = path_merged_data + "final-merged-dfs.txt"
    params:
        ds_param_dict_list = config["ds_param_dict_list"],
        build_datatable = config["build_datatable"]
    script:
        "scripts/integrate_datasets.R"

# Create "delta datasets"; differences between samples from different
# timepoints, within individual
rule make_delta_datasets:
    input:
        in_file = path_merged_data + "final-merged-dfs.txt"
    output:
        out_file = config["out"] + "04-SELECTED-FEATURES-{reference}/"
                                   "deltas-{reference}.txt"
    params:
        reference_time = "{reference}",
        absolute_values = "no",
        build_datatable = config["build_datatable"], # TODO have this be optional default true
        distance_matrices = config["distance_matrices"]
    script:
        "scripts/create_deltas.R"


# rule visualize_datasets:
#     input:
#         in_file = config["out"] + "04-SELECTED-FEATURES-{reference}/deltas-{reference}.txt"
#     output:
#         out_file = config["out"] + "04-SELECTED-FEATURES-{reference}/vizualizer-{reference}.html"
#     script:
#         "scripts/viz-datatable.R"


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
rule random_forest_deltas:
    input:
        in_file = config["out"] +
                  "04-SELECTED-FEATURES-{reference}/deltas-{reference}.txt"
    output:
        out_file = config["out"] +
                   "04-SELECTED-FEATURES-{reference}/"
                   "{mixed}-RF-{deltas}-{reference}.pdf"
    params:
        random_forest_type = "{mixed}",  # mixed or fixed
        random_effect = config["random_effect"],
        sample_id = config["sample_id"],
        response_var = config["response_var"],
        delta_flag = "{deltas}",  # raw or deltas
        iterations = config["iterations"],  # 20 is suggested, 10 for testing
    conda: "conda_envs/merf.yaml"
    script:
        "scripts/random_forest.py"


rule random_forest_original:
    input:
        in_file = path_merged_data + "final-merged-dfs.txt"
    output:
        out_file = path_rf_original + "MERF-original.pdf"
    params:
        random_forest_type = "mixed",
        random_effect = config["random_effect"],
        sample_id = config["sample_id"],
        response_var = config["response_var"],
        delta_flag = "raw",
        iterations = config["iterations"],
    conda: "conda_envs/merf.yaml"
    script:
        "scripts/random_forest.py"


# TODO program if statements here and make report reflect missing data
# TODO consider if statements with all rule
# TODO idea make four new render reports
# TODO consider making some of these parameters
rule run_post_hoc_stats:
    input:
        deltas_first = path_rf_first + "deltas-first.txt",
        deltas_previous = path_rf_previous + "deltas-previous.txt",
        deltas_pairwise = path_rf_pairwise + "deltas-pairwise.txt",
        df_original = path_merged_data + "final-merged-dfs.txt",
        fixed_effects_first = path_rf_first + "mixed-RF-deltas-first-boruta-important.txt",
        fixed_effects_previous = path_rf_previous + "mixed-RF-deltas-previous-boruta-important.txt",
        fixed_effects_pairwise = path_rf_pairwise + "mixed-RF-deltas-pairwise-boruta-important.txt",
        fixed_effects_original = path_rf_original + "MERF-original-boruta-important.txt"
    output:
        out_file = path_post_hoc + "post-hoc-analysis.html"
    params:
        response_var = "HDL"
    script:
        "scripts/post_hoc_tests.R"


rule render_report:
    input:
        original = rules.random_forest_original.output.out_file,
        original_shap= path_rf_original + "MERF-original-accepted-SHAP-"
                                          "summary-beeswarm.svg",
        original_boruta = path_rf_original + "MERF-original-boruta-"
                        "accepted-features.svg",
        original_log = path_rf_original + "MERF-original-log.txt",

        first = path_rf_first + "mixed-RF-deltas-first.pdf",
        first_shap= path_rf_first + "mixed-RF-deltas-first-accepted-SHAP-"
                                    "summary-beeswarm.svg",
        first_boruta= path_rf_first + "mixed-RF-"
                      "deltas-first-boruta-accepted-features.svg",
        first_log = path_rf_first + "mixed-RF-deltas-first-log.txt",

        previous = path_rf_previous + "mixed-RF-deltas-previous.pdf",
        previous_shap= path_rf_previous + "mixed-RF-deltas-previous-"
                                          "accepted-SHAP-summary-beeswarm.svg",
        previous_boruta = path_rf_previous + "mixed-RF-deltas-previous-"
                                             "boruta-accepted-features.svg",
        previous_log= path_rf_previous + "mixed-RF-deltas-previous-log.txt",

        pairwise = path_rf_pairwise + "mixed-RF-deltas-pairwise.pdf",
        pairwise_shap= path_rf_pairwise + "mixed-RF-deltas-pairwise-"
                                          "accepted-SHAP-summary-beeswarm.svg",
        pairwise_boruta = path_rf_pairwise + "mixed-RF-deltas-pairwise-"
                                             "boruta-accepted-features.svg",
        pairwise_log = path_rf_pairwise + "mixed-RF-deltas-pairwise-log.txt",
        post_hoc = path_post_hoc + "post-hoc-analysis.html"
    output:
        md_doc=config["report_name"]
    script:
        "scripts/report.Rmd"

