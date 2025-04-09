
import os

import shap
import pandas as pd
import numpy as np
import sys
import re

import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from matplotlib.backends.backend_pdf import PdfPages

from merf.merf import MERF
from boruta_shap_min.borutashap import BorutaShap
from boruta_shap_plot import _boruta_shap_plot

from first_checks import first_checks


# from config
response_var = snakemake.params["response_var"]

# from rule
dataset = snakemake.params["dataset"]

shap_values = snakemake.input["shap_values"]
features_for_shap = snakemake.input["features_for_shap"]
important_features = snakemake.input["important_features"]

out_file = snakemake.output["out_file"]

# set graphics dimensions
_width = 10
_height = 6

join_flag = True

ds_out_path = snakemake.config["out"] + snakemake.config["path_" + dataset]

plot_list_shap = []
plot_file_name_list_shap = []

plot_list_dependence = []
plot_file_name_list_dependence = []


class Logger(object):
    def __init__(self):
        self.log = open(ds_out_path + dataset + "-log.txt", "a")

    def write(self, message):
        self.log.write(message)

    def flush(self):
        pass


def _build_result_pdf(out_file, out_file_prefix, plot_list,
                      plot_file_name_list):
    with PdfPages(out_file) as pdf:
        # TODO add only important plot list (don't save svg or pdf)
        for i, plt in enumerate(plot_list):
            file_name = out_file_prefix + plot_file_name_list[i]

            plt.savefig(file_name + ".svg", format='svg')

            pdf.savefig(plt)


def graphic_handling(plt, plot_file_name, p_title, width,
                     height, tick_size, title_size, plot_group,
                     shap_statement):
    if plot_group == "shap" or plot_group == "dependence":
        plot_file_name_list_shap.append(plot_file_name)
    if plot_group == "dependence":
        plot_file_name_list_dependence.append(plot_file_name)

    top_adjustment = 0.80 if plot_group == "shap" else 0.95

    plt.title(p_title, fontdict={"fontsize": title_size})
    plt.tick_params(axis="both", labelsize=tick_size)
    plt.annotate(shap_statement, xy=(top_adjustment, 1.08),
                 xycoords='axes fraction', fontsize=6, ha='center',
                 va='center')
    p = plt.gcf()
    p.set_size_inches(width, height)
    plt.close("all")

    if plot_group == "shap" or plot_group == "dependence":
        plot_list_shap.append(p)
    if plot_group == "dependence":
        plot_list_dependence.append(p)


def list_split(input_list, group_size):
    max_items = group_size
    return [input_list[i:i+max_items]
            for i in range(0, len(input_list), max_items)]


def my_shap_summary_plots(feature_index_list, shap_values, x, list_num,
                          group_size, shap_statement):
    plt_s = shap_values[:, feature_index_list]
    plt_x = x.iloc[:, feature_index_list]
    ax_sz = 8
    title_sz = 8
    tick_sz = 6

    # beeswarm
    plt.figure(figsize=(_width, _height))
    d = str.capitalize(dataset)
    t = f"Features related to {response_var} ({d} dataset)"
    shap.summary_plot(shap_values=plt_s, features=plt_x, show=False,
                      sort=True, max_display=group_size, color_bar=True,
                      color_bar_label='Selected feature value',
                      plot_size=(0.8*_width, 0.7*_height))
    plt.ylabel("Selected Features\n", fontsize=ax_sz)
    plt.xlabel("SHAP value (impact on response)\nLeft of zero is " +
               "negative impact and right of zero is positive impact",
               fontsize=ax_sz)
    graphic_handling(plt, "-accepted-SHAP-summary-beeswarm-" +
                     str(list_num), t, _width, _height, tick_sz,
                     title_sz, "shap", shap_statement)


def shap_plots(shap_values, x, boruta_accepted, feature_importance,
               shap_statement):

    # get indices from Boruta for accepted features
    index_list_accepted = feature_importance["feature_index"].tolist()

    # Number of important features can vary, so split up graphs so
    # graphs aren't cluttered (note scales will vary)
    group_size = 10
    lists_to_graph = list_split(index_list_accepted, group_size=group_size)

    # graph each list of important features after
    for i in range(0, len(lists_to_graph)):
        my_shap_summary_plots(lists_to_graph[i], shap_values, x, i,
                              group_size, shap_statement)

    col_count = 1

    for col in boruta_accepted:
        # Dependence plots
        plt.figure(figsize=(_width, _height))
        shap.dependence_plot(col, shap_values, x, show=False, x_jitter=0.03,
                             interaction_index="auto")
        plt.ylabel("SHAP value for\n" + str(col), fontsize=9)
        plt.xlabel(str(col), fontsize=9)
        graphic_handling(plt, "-SHAP-dependence-plot-interaction-" +
                         str(col_count),
                         "Dependence plot showing feature impact on\n " +
                         response_var + " (" + str.capitalize(dataset) +
                         " dataset)", _width, _height, 7, 8,
                         "dependence", shap_statement)
        col_count += 1


def _cleanup_files():
    # some files were used just to compile pdf, so remove them
    out_path = snakemake.config["out"] + snakemake.config["path_" + dataset]
    files_in_folder = os.listdir(out_path)
    for file_name in files_in_folder:
        for n in ["SHAP-dependence-plot", "features.svg"]:
            if n in file_name:
                fp = os.path.join(out_path, file_name)
                os.remove(fp)


def main(shap_values, features_for_shap, important_features, out_file):

    sys.stdout = Logger()

    if snakemake.config["analyze_" + dataset] == "no":
        with open(out_file, "w") as file:
            file.write("")
        return

    # if response var is categorical, needs to be factorized and
    # mapped to numeric values to be used with RF regressor
    shap_statement = ""
    df = pd.read_csv(ds_out_path + dataset + ".txt", sep="\t")
    if df[response_var].dtype == "object":
        df = df.sort_values(by=response_var)
        df[response_var], u_maps = pd.factorize(df[response_var])

        print("\nNumeric mapping was created for categoric response " +
              "variable, " + str(response_var) +
              ", (sorted alphabetically and factorized):")
        shap_statement = str(response_var) + " (response) mapping:\n"
        for i, category in enumerate(u_maps):
            print(f"{category} maps to {i}")
            shap_statement += f"{category} maps to {i}\n"
        print("\n*** Percent variation explained is not optimal for " +
              "categorical response variables. \n*** Use caution with " +
              "interpretation.\n")

    shap_values_df = pd.read_csv(shap_values, sep="\t")
    shap_values_arr = shap_values_df.values

    important_features = pd.read_csv(important_features, sep="\t")
    important_features_names = important_features['important_features']

    features_for_shap = pd.read_csv(features_for_shap, sep="\t")

    shap_plots(shap_values_arr, features_for_shap, important_features_names,
               important_features, shap_statement)

    plt.close("all")
    # PDF; Primary figures (SHAP)
    _build_result_pdf(out_file, ds_out_path + dataset, plot_list_shap,
                      plot_file_name_list_shap)
    _build_result_pdf(ds_out_path + dataset + "-dependence.pdf",
                      ds_out_path + dataset, plot_list_dependence,
                      plot_file_name_list_dependence)

    _cleanup_files()


main(shap_values=shap_values, features_for_shap=features_for_shap,
     important_features=important_features, out_file=out_file)
