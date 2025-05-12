
import os
import pandas as pd


def _first_checks(config):
    required = ["analyst", "random_effect", "sample_id", "timepoint",
                "include_time", "out", "response_var", "analyze_first",
                "analyze_previous", "analyze_pairwise", "analysis_notes",
                "iterations", "n_estimators", "max_features",
                "absolute_values", "include_reference_values",
                "borutashap_trials", "borutashap_threshold", "borutashap_p",
                "enc_percent_threshold", "distance_matrices", "input_datasets",
                "path_dim_pca", "path_dim_scnic", "path_original",
                "path_first", "path_previous", "path_pairwise", "version",
                "file_path_list", "ds_param_dict_list", "pca_ds_list",
                "scnic_ds_list"]

    if sorted(required) == sorted(config.keys()):
        with open(f"{config['out']}check.done", 'w') as file:
            file.write("")
    else:
        print("Missing arguments in config file")


def write_files(path, dataset, config, summary_table_items):

    repeat_na = len(summary_table_items) - 1
    blank = {
        "Data": summary_table_items,
        str.capitalize(dataset):  ["Not performed"] + (["NA"] * repeat_na)
    }
    blank_df = pd.DataFrame(blank)

    boruta_important_df = pd.DataFrame({"important_features":
                                        ["no_selected_features"],
                                        "decoded_features": ["NA"],
                                        "feature_importance_vals": [-100]})
    boruta_important_df = pd.DataFrame(boruta_important_df)

    out_files1 = ["", "", "", "", ""]
    # create blank files for when Original, First, Previous or Pairwise
    # analyses are not performed
    extensions1 = ["-boruta.pdf", "-log.txt", "-input-model-df.txt",
                   "-SHAP-values-df.txt", "-features-for-shap.txt"]


    for ext, content in zip(extensions1, out_files1):
        out_file = f"{config['out']}{config[path]}{dataset}{ext}"
        with open(out_file, "w") as file:
            file.write(content)

    # files with dataframes
    out_files2 = [boruta_important_df, blank_df]
    extensions2 = ["-boruta-important.txt", "-log-df.txt"]

    for ext, content in zip(extensions2, out_files2):
        out_file = f"{config['out']}{config[path]}{dataset}{ext}"
        content.to_csv(out_file, index=False, sep="\t")


def main(config, dataset, summary_table_items):
    config = config

    dataset_path = "path_" + dataset

    dir_name = config["out"] + config[dataset_path]
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    analyze_current_dataset = config["analyze_" + dataset]

    if analyze_current_dataset == "no":
        write_files(dataset_path, dataset, config, summary_table_items)

    # _first_checks(config) # use to check for required items


def first_checks(config, dataset, summary_table_items):
    main(config, dataset, summary_table_items)
