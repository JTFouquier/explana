
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

# From config
borutashap_trials = int(snakemake.config["borutashap_trials"])
borutashap_threshold = int(snakemake.config["borutashap_threshold"])
borutashap_p = float(snakemake.config["borutashap_p"])
enc_percent_threshold = float(snakemake.config["enc_percent_threshold"])

# From rule
random_effect = snakemake.params["random_effect"]
sample_id = snakemake.params["sample_id"]
response_var = snakemake.params["response_var"]
include_time = str.lower(snakemake.params["include_time"])

n_estimators = int(snakemake.params["n_estimators"])
max_features = float(snakemake.params["max_features"])
max_depth = 7  # depth recommended by BorutaSHAP; kept constant

iterations = int(snakemake.params["iterations"])

df_input = snakemake.input["in_file"]
out_file = snakemake.output["out_file"]
dataset = snakemake.params["dataset"]

# set graphics dimensions
_width = 10
_height = 6

join_flag = True

ds_out_path = snakemake.config["out"] + snakemake.config["path_" + dataset]

plot_list_shap = []
plot_file_name_list_shap = []

plot_list_dependence = []
plot_file_name_list_dependence = []

plot_list_boruta = []
plot_file_name_list_boruta = []

summary_table_items = ["Model Type (Pass/Fail)", "% Variance Explained",
                       "N Trees", "Feature fraction/split", "Max Depth",
                       "MERF Iters.", "BorutaSHAP Trials",
                       "BorutaSHAP Threshold", "P-value", "N Study IDs",
                       "N Samples", "Input Features", "Accepted Features",
                       "Tentative Features", "Rejected Features"]


# saves to log file; does not print to screen
class Logger(object):
    def __init__(self):
        self.log = open(ds_out_path + dataset + "-log.txt", "a")

    def write(self, message):
        self.log.write(message)

    def flush(self):
        pass


def _build_result_pdf(out_file, out_file_prefix, plot_list, plot_file_name_list):
    with PdfPages(out_file) as pdf:
        # TODO add only important plot list (don't save svg or pdf)
        for i, plt in enumerate(plot_list):
            file_name = out_file_prefix + plot_file_name_list[i]

            plt.savefig(file_name + ".svg", format='svg')

            pdf.savefig(plt)


# placeholder dataframe for creating figure
def no_features_selected():
    _df = pd.DataFrame({"important_features": ["no_selected_features"],
                        "decoded_features": ["NA"],
                        "feature_importance_vals": [-100]})
    _df.to_csv(ds_out_path + dataset + "-boruta-important.txt", index=False,
               sep="\t")


def fail_analysis(out_file):
    with open(out_file, "w") as file:
        file.write("Failed Analysis")
    no_features_selected()
    return


def dummy_var_creation(categorical_list, df, join_flag):
    print("\nBinary encoded columns created for categorical input variables:")
    for cat_var in categorical_list:
        subset_df = df[[cat_var]]
        encoded_name = "ENC_" + cat_var + "_is"
        dum_df = pd.get_dummies(subset_df, columns=[cat_var],
                                prefix=[encoded_name])
        print("-> " + cat_var + ": " + str(dum_df.columns.to_list()))
        # Keep adding encoded vars to df before returning
        if join_flag:
            df = df.join(dum_df)
    if join_flag:
        return df

    else:
        return dum_df


def setup_df_do_encoding(df, random_effect, response_var, join_flag,
                         sample_id):
    shap_statement = ""  # Not needed unless categorical response
    numeric_column_list = df._get_numeric_data().columns.values
    categoric_columns = [i for i in list(df.columns) if
                         i not in numeric_column_list]

    if include_time == "yes":
        do_not_encode = [random_effect, sample_id, response_var]
    if include_time == "no":
        do_not_encode = [random_effect, sample_id, response_var,
                         "timepoint_explana"]

    study_n = str(df[random_effect].nunique())

    # encode variables to improve looking at specific values
    vars_to_encode = [i for i in categoric_columns if i not in do_not_encode]
    df = dummy_var_creation(vars_to_encode, df, join_flag)

    # if response var is categoric convert to integers
    if df[response_var].dtype == "object":

        df = df.sort_values(by=response_var)

        # remove vars before encoding, random effect, response, sample ID
        df[response_var], u_maps = pd.factorize(df[response_var])

        print("\nNumeric mapping was created for categoric response" +
              "variable, " + str(response_var) +
              ", (sorted alphabetically and factorized):")
        shap_statement = str(response_var) + " (response) mapping:\n"
        for i, category in enumerate(u_maps):
            print(f"{category} maps to {i}")
            shap_statement += f"{category} maps to {i}\n"
        print(shap_statement)
        print("\n*** Percent variation explained is not optimal for " +
              "categorical response variables. \n*** Use caution with " +
              "interpretation.\n")

    # always drop original because timepoint_explana is ranked
    if include_time == "yes":
        df_encoded = df.drop(do_not_encode + vars_to_encode +
                             ["timepoint_original"], axis=1)

    if include_time == "no":
        df_encoded = df.drop(do_not_encode + vars_to_encode +
                             ["timepoint_original", "timepoint_explana"],
                             axis=1)

    def remove_low_percentage_encoded_vars(df, percent):
        # TODO make function for removing columns
        # yes/1 counts per column; drop columns in
        non_encoded_cols = [col for col in df.columns if not "ENC_" in col]
        df_categoric = df.drop(non_encoded_cols, axis=1)

        value_counts = df_categoric.eq(1).sum()
        target_value = round(len(df_categoric)*(percent/100))

        # list of column names to drop from df
        columns_to_drop = list(value_counts[value_counts < target_value].index)
        df = df.drop(columns_to_drop, axis=1)
        feature_list = df.columns.to_list()

        return (feature_list)

    if enc_percent_threshold >= 0:
        feature_list = \
            remove_low_percentage_encoded_vars(df=df_encoded,
                                               percent=enc_percent_threshold)
    else:
        feature_list = df.columns.to_list()

    x = df[feature_list]  # Raw analysis and ALL (new, remove ref)
    z = np.ones((len(x), 1))
    c = df[random_effect]
    y = df[response_var]

    x.to_csv(ds_out_path + dataset + "-input-model-df.txt", index=False,
             sep="\t")
    df_features_in = pd.DataFrame(x.columns, columns=["input_features"])

    return x, z, c, y, df_features_in, study_n, shap_statement


def train_rf_model(x, y, rf_regressor, step):
    forest = rf_regressor.fit(x, y)
    oob = str(round(forest.oob_score_*100, 1))  # percent variation
    plt.close("all")
    return forest, oob


def graphic_handling(plt, plot_file_name, p_title, width,
                     height, tick_size, title_size, plot_group, shap_statement):
    if plot_group == "shap" or plot_group == "dependence":
        plot_file_name_list_shap.append(plot_file_name)
    if plot_group == "boruta":
        plot_file_name_list_boruta.append(plot_file_name)
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
    if plot_group == "boruta":
        plot_list_boruta.append(p)
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
    index_list_accepted = list(feature_importance[feature_importance[
        "col_name"].isin(boruta_accepted)].index)

    # Number of important features can vary, so split up graphs so
    # graphs aren't cluttered (note scales will vary)
    group_size = 10
    lists_to_graph = list_split(index_list_accepted, group_size=group_size)

    shap_values_df = pd.DataFrame(shap_values, columns=list(x.columns.values))
    shap_values_df.to_csv(ds_out_path + dataset + "-SHAP-values-df.txt",
                          index=False, sep="\t")
    # graph each list of important features after
    for i in range(0, len(lists_to_graph)):
        my_shap_summary_plots(lists_to_graph[i], shap_values, x, i, group_size,
                              shap_statement)

    col_count = 1
    # TODO save the dataframe and get vals 1
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


def run_shap(trained_forest_model, x):
    shap_explainer = \
        shap.TreeExplainer(trained_forest_model,
                           feature_perturbation="tree_path_dependent")
    shap_values = shap_explainer.shap_values(x)

    feature_names = list(x.columns.values)
    vals = np.abs(shap_values).mean(0)
    shap_imp = \
        pd.DataFrame(list(zip(feature_names, vals)), columns=["col_name",
                     "feature_importance_vals"])
    shap_imp.sort_values(by=["feature_importance_vals"],
                         ascending=False, inplace=True)
    shap_imp["feature_importance_vals"] = \
        shap_imp["feature_importance_vals"].apply(lambda x: round(x, 4))

    return shap_imp, shap_values


def run_boruta_shap(forest, x, y, shap_statement):
    boruta_height = 10
    feature_selector = \
        BorutaShap(model=forest, importance_measure="shap",
                   classification=False, percentile=borutashap_threshold,
                   pvalue=borutashap_p)
    feature_selector.fit(X=x, y=y, n_trials=borutashap_trials,
                         train_or_test="train", sample=False, verbose=False)
    for feature_group in ["all", "accepted", "tentative", "rejected"]:
        boruta_width = _boruta_shap_plot(feature_selector,
                                         which_features=feature_group)
        graphic_handling(plt, "-boruta-" + feature_group + "-features",
                         "Boruta " + feature_group + " features",
                         boruta_width, boruta_height, 5, 9, "boruta",
                         shap_statement)
    return feature_selector


def train_merf_model(x, z, c, y, model, step):
    mrf = MERF(max_iterations=iterations, fixed_effects_model=model)
    mrf.fit(x, z, c, y)
    forest = mrf.trained_fe_model
    oob = str(round(forest.oob_score_*100, 1))  # percent variation
    return forest, oob


def summary_stats_table(x, df_boruta):
    x_subset = x[df_boruta["important_features"]]
    feat_names = x_subset.columns
    x_subset = x_subset.describe(include='all').transpose()
    x_subset = x_subset.round(3)

    x_subset["important_features"] = feat_names
    merged_df = pd.merge(x_subset, df_boruta, on='important_features',
                         how='left')
    del merged_df["decoded_features"]
    del merged_df["count"]
    first = ["important_features", "feature_importance_vals"]
    my_order = first + [col for col in merged_df if col not in first]
    merged_df = merged_df[my_order]

    merged_df.to_csv(ds_out_path + dataset + "-feature-stats.txt", index=False,
                     sep="\t")


def _cleanup_files():
    # some files were used just to compile pdf, so remove them
    out_path = snakemake.config["out"] + snakemake.config["path_" + dataset]
    files_in_folder = os.listdir(out_path)
    for file_name in files_in_folder:
        for n in ["SHAP-dependence-plot", "features.svg"]:
            if n in file_name:
                fp = os.path.join(out_path, file_name)
                os.remove(fp)


def main(df_input, out_file, random_effect, sample_id, response_var,
         join_flag):

    # make a log for printed information and a dataframe
    sys.stdout = Logger()

    first_checks(snakemake.config, dataset, summary_table_items)

    # TODO this should not need to happen... should not have to initiate
    # the environment just to find out the user didn't want to build model
    if snakemake.config["analyze_" + dataset] == "no":
        return

    # for model summary table
    model_evaluation = "NA"
    var_explained = "NA"
    study_n = "NA"
    sample_n = "NA"
    input_features = "NA"
    accepted_features = "NA"
    tentative_features = "NA"
    rejected_features = "NA"
    merf_iters = "NA"

    summary = {
        "Data": summary_table_items,
        str.capitalize(dataset): ["NA"] * len(summary_table_items)
    }
    summary = pd.DataFrame(summary)

    df = pd.read_csv(df_input, sep="\t", na_filter=False)

    # if the input dataframe is empty; rule terminates
    rf_type = ""
    if "Failed Analysis" in df.columns:
        model_evaluation = "Not performed"
        fail_analysis(out_file)
        fail_analysis(ds_out_path + dataset + "-boruta.pdf")
        fail_analysis(ds_out_path + dataset + "-dependence.pdf")

    else:
        study_n = str(df[random_effect].nunique())
        sample_n = str(df[sample_id].nunique())

        if len(np.unique(df[random_effect])) != len(df[random_effect]):
            rf_type = "MERF"
            merf_iters = str(int(iterations))
        else:
            rf_type = "RF"
            merf_iters = "NA"

        x, z, c, y, df_features_in, study_n, shap_statement = \
            setup_df_do_encoding(df=df, random_effect=random_effect,
                                 response_var=response_var,
                                 join_flag=join_flag, sample_id=sample_id)

        # Set up RF regressor for both fixed or mixed models
        def instantiate_rf():
            rf_regressor = \
                RandomForestRegressor(n_jobs=-1, n_estimators=n_estimators,
                                      oob_score=True,
                                      max_features=max_features,
                                      max_depth=max_depth)
            return (rf_regressor)

        # Run either mixed effects or fixed effects RF models
        if rf_type == "MERF":
            forest, oob = train_merf_model(x=x, z=z, c=c, y=y,
                                           model=instantiate_rf(),
                                           step="full_forest")
            plt.close("all")

        elif rf_type == "RF":
            forest, oob = train_rf_model(x, y, instantiate_rf(),
                                         step="full_forest")

        # if first time running RF is < 5.0 percent variation, end analysis
        var_explained = str(oob) + "%"
        if float(oob) < 5.0:
            model_evaluation = rf_type + " FAIL (low variance)"
            fail_analysis(out_file)
            fail_analysis(ds_out_path + dataset + "-boruta.pdf")
            fail_analysis(ds_out_path + dataset + "-dependence.pdf")

        else:
            model_evaluation = rf_type + " PASS"
            # BorutaSHAP to select features
            boruta_explainer = run_boruta_shap(forest, x, y, shap_statement)
            # refit forest after excluding rejected features prior to
            # visualizing data. Essentially if data is visualized in SHAP plots
            # with all the rejected features, they're not effective
            # for understanding feature impact on response
            # z, c, and y are the same as before

            rejected_cols = boruta_explainer.rejected

            if "timepoint_explana" in boruta_explainer.rejected:
                rejected_cols.remove("timepoint_explana")

            x = x.drop(rejected_cols, axis=1)

            if rf_type == "MERF":
                forest, oob_rerun = train_merf_model(x=x, z=z, c=c, y=y,
                                                     model=instantiate_rf(),
                                                     step="shap_rerun")
            elif rf_type == "RF":
                forest, oob_rerun = train_rf_model(x, y, instantiate_rf(),
                                                   step="shap_rerun")

            var_explained += " (" + str(oob_rerun) + "%)"
            shap_imp, shap_values = run_shap(forest, x)
            shap_plots(shap_values, x, boruta_explainer.accepted,
                       shap_imp, shap_statement)
            # get encoded vals prefix "ENC_Location_is_1_3" decoded = Location
            decoded_top_features_list = []
            shap_imp_score_list = []
            # get original colname if feature was encoded; else use var name
            for feature_name in boruta_explainer.accepted:
                try:
                    r = re.search("ENC_(.*)_is_", feature_name)
                    decoded = r.group(1)
                except AttributeError:
                    decoded = feature_name

                decoded_top_features_list.append("".join(decoded))
                shap_imp_score = shap_imp.loc[shap_imp["col_name"] ==
                                              feature_name,
                                              "feature_importance_vals"]
                shap_imp_score_list.append(shap_imp_score.iloc[0])

            df_boruta = pd.DataFrame({
                "important_features": boruta_explainer.accepted,
                "decoded_features": decoded_top_features_list,
                "feature_importance_vals": shap_imp_score_list})
            df_boruta = df_boruta.sort_values(by=["feature_importance_vals"],
                                              ascending=False)

            df_boruta.to_csv(ds_out_path + dataset +
                             "-boruta-all-features.txt",
                             index=False, sep="\t")

            # check if input feature was selected or not
            def check_inputs(row):
                if pd.isnull(row["important_features"]):
                    return "no"
                else:
                    return "yes"

            df_features_in = \
                df_features_in.merge(df_boruta["important_features"],
                                     how="left", left_on="input_features",
                                     right_on="important_features")
            df_features_in["was_selected"] = df_features_in.apply(check_inputs,
                                                                  axis=1)
            df_features_in = df_features_in.drop(["important_features"],
                                                 axis=1)
            df_features_in.to_csv(ds_out_path + dataset +
                                  "-input-features.txt", sep="\t",
                                  index=False)

            # data for log data frame for report
            input_features = str(len(boruta_explainer.all_columns))
            accepted_features = str(len(boruta_explainer.accepted))
            tentative_features = str(len(boruta_explainer.tentative))
            rejected_features = str(len(boruta_explainer.rejected))

            if df_boruta.empty:
                # Make model summary table for final report
                model_evaluation = "FAIL: " + model_evaluation + "; Boruta FAIL"
                summary[str.capitalize(dataset)] = [model_evaluation,
                                                    var_explained,
                                                    n_estimators, max_features,
                                                    max_depth, merf_iters,
                                                    borutashap_trials,
                                                    borutashap_threshold,
                                                    borutashap_p, study_n,
                                                    sample_n, input_features,
                                                    accepted_features,
                                                    tentative_features,
                                                    rejected_features]

                summary.to_csv(ds_out_path + dataset + "-log-df.txt",
                               index=False, sep="\t")
                fail_analysis(out_file)
                plt.close("all")
                # PDF Boruta figs; helpful even if none selected by Boruta
                _build_result_pdf(ds_out_path + dataset + "-boruta.pdf",
                                  ds_out_path + dataset, plot_list_boruta,
                                  plot_file_name_list_boruta)
                _cleanup_files()
                print('No features selected', file=sys.stderr)
                sys.exit(0)
            else:
                model_evaluation = model_evaluation + "; Boruta PASS"
                df_boruta.to_csv(ds_out_path + dataset +
                                 "-boruta-important.txt",
                                 index=False, sep="\t")

                summary_stats_table(x, df_boruta)

            plt.close("all")
            # PDF; Primary figures (SHAP)
            _build_result_pdf(out_file, ds_out_path + dataset, plot_list_shap,
                              plot_file_name_list_shap)
            # PDF; Boruta figures
            _build_result_pdf(ds_out_path + dataset + "-boruta.pdf",
                              ds_out_path + dataset, plot_list_boruta,
                              plot_file_name_list_boruta)
            # PDF; Dependence plots
            _build_result_pdf(ds_out_path + dataset + "-dependence.pdf",
                              ds_out_path + dataset, plot_list_dependence,
                              plot_file_name_list_dependence)

        # Make model summary table for final report
        summary[str.capitalize(dataset)] = [model_evaluation, var_explained,
                                            n_estimators, max_features,
                                            max_depth, merf_iters,
                                            borutashap_trials,
                                            borutashap_threshold, borutashap_p,
                                            study_n, sample_n,
                                            input_features, accepted_features,
                                            tentative_features,
                                            rejected_features]

        summary.to_csv(ds_out_path + dataset + "-log-df.txt", index=False,
                       sep="\t")

    _cleanup_files()


main(df_input=df_input, out_file=out_file, random_effect=random_effect,
     sample_id=sample_id, response_var=response_var, join_flag=join_flag)
