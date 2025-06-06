import json
import os
import sys
import shutil

from snakemake.common.configfile import load_configfile

# add vars that users likely don't need
# if they need them, they'll be loaded next
config["delta_df_mod"] = "delta_df <- delta_df"

if not os.path.exists(config["out"]):
    os.makedirs(config["out"])

index = sys.argv.index("--configfile")
config_file = sys.argv[index + 1]
out_dir = config["out"]

# make a copy of the config file for repository
new_filename = "analysis_config.yaml"
shutil.copy(config_file, os.path.join(out_dir, new_filename))

# After loading config; don't want users modifying paths
config["path_dim_pca"] = "DIM-PCA/"
config["path_dim_scnic"] = "DIM-SCNIC/"
config["path_dim_preprocess"] = "DIM-preprocess/"

# can simplifity this code
config["path_original"] = "SELECTED-FEATURES-original/"
config["path_first"] = "SELECTED-FEATURES-first/"
config["path_previous"] = "SELECTED-FEATURES-previous/"
config["path_pairwise"] = "SELECTED-FEATURES-pairwise/"

config["version"] = "2025.05.09"

summary_table_items = ["Model Type (Pass/Fail)", "% Variance Explained",
                       "N Trees", "Feature fraction/split", "Max Depth",
                       "MERF Iters.", "BorutaSHAP Trials",
                       "BorutaSHAP Threshold", "P-value", "N Study IDs",
                       "N Samples", "Input Features", "Accepted Features",
                       "Tentative Features", "Rejected Features"]

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
path_dim_preprocess = out_dir + config["path_dim_preprocess"]


def dim_reduction():
    input_datasets = config["input_datasets"]
    file_path_list = []
    ds_param_dict_list = []
    pca_ds_list = []
    scnic_ds_list = []
    preprocess_ds_list = []

    for ds_name in input_datasets:
        dim_method = input_datasets[ds_name]["dim_method"]
        ds_param_dict_list.append(input_datasets[ds_name]["df_mod"])

        dim_method = dim_method.lower()
        if dim_method == "pca":
            f_out = path_dim_pca + ds_name + "/PCA-dim-reduction-" \
                                             "for-workflow.txt"
            pca_ds_list.append(ds_name)

        elif dim_method == "scnic":
            f_out = path_dim_scnic + ds_name + \
                    "/SCNIC-dim-reduction-for-workflow.txt"
            scnic_ds_list.append(ds_name)

        elif dim_method == "preprocess" or dim_method == "":
            f_out = path_dim_preprocess + ds_name + \
            "/preprocess-for-workflow.txt"
            preprocess_ds_list.append(ds_name)

        elif dim_method == "none":
            f_out = input_datasets[ds_name]["file_path"]

        file_path_list.append(f_out)

    return file_path_list, ds_param_dict_list, \
           pca_ds_list, \
           scnic_ds_list, preprocess_ds_list, input_datasets

config["file_path_list"], config["ds_param_dict_list"], \
config["pca_ds_list"], \
config["scnic_ds_list"], config["preprocess_ds_list"], input_datasets = dim_reduction()
config["input_datasets"] = input_datasets


def get_preprocess_df_mod(wildcards):
    df_mod = input_datasets[wildcards.preprocess_ds_name]["df_mod"]
    return df_mod

def get_preprocess_in_file(wildcards):
    in_file = input_datasets[wildcards.preprocess_ds_name]["file_path"]
    return in_file

def get_preprocess_method(wildcards):
    method = input_datasets[wildcards.preprocess_ds_name]["dim_param_dict"]["method"]
    return method

# def get_pca_df_mod(wildcards):
#     df_mod = input_datasets[wildcards.pca_ds_name]["df_mod"]
#     return df_mod

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


rule preprocess:
    input:
        in_file = get_preprocess_in_file,
    output:
        out_file = path_dim_preprocess + "{preprocess_ds_name}" +
        "/preprocess-for-workflow.txt"
    params:
        dataset_name = expand("{preprocess_ds_name}",
            preprocess_ds_name=config["preprocess_ds_list"]),
        method = get_preprocess_method,
        df_mod = get_preprocess_df_mod
    conda: "conda_envs/r_env.yaml",
    script:
        "scripts/preprocess.R"


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
        out_boruta = protected(out_dir + \
        "SELECTED-FEATURES-{reference}/{reference}-boruta-important.txt"),
        model_df = protected(out_dir + "SELECTED-FEATURES-{reference}/{reference}-input-model-df.txt"),
        shap_values = protected(out_dir + "SELECTED-FEATURES-{reference}/{reference}-SHAP-values-df.txt"),
        features_for_shap = protected(out_dir + "SELECTED-FEATURES-{reference}/{reference}-features-for-shap.txt"),
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


rule shap_plots:
    input:
        shap_values = out_dir + "SELECTED-FEATURES-{reference}/{reference}-SHAP-values-df.txt",
        features_for_shap = out_dir + "SELECTED-FEATURES-{reference}/{reference}-features-for-shap.txt",
        important_features = out_dir + "SELECTED-FEATURES-{reference}/{reference}-boruta-important.txt"
    output:
        out_file = out_dir + "SELECTED-FEATURES-{reference}/{reference}.pdf",
    params:
        dataset = "{reference}",
        response_var = config["response_var"],
    conda: "conda_envs/merf.yaml"
    script:
        "scripts/shap_plots.py"


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