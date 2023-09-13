

# TODO organize imports; these are a mess
import os

import shap
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import re

from merf.merf import MERF
from BorutaShap import BorutaShap
# from boruta_shap_mine import BorutaShap
from sklearn.ensemble import RandomForestRegressor

# TODO switch to OneHotEncoder
from utils import _build_result_pdf
from boruta_shap_plot import _boruta_shap_plot

"""
TODO list 
Add baseline information and perform tests

Simple fixes
- add use reference for deltas with categorical 
- fix hardcoding & add sample ID for all
- consider shadowed scope vars
- fix rounding for feature importance output file
- remove leading 'X' from feature names

Features/Considerations for development
- consider adding time as both categorical and continuous variable
    and allow flexibility
- consider what random effects make sense -- allow user to choose? 
- longitudinal vs time series
- gradient descent (ML) vs inverse matrix computation formula (stats) 
    (for one-hot encoding)

# reading
# https://slundberg.github.io/shap/notebooks/plots/dependence_plot.html
# https://shap.readthedocs.io/en/latest/example_notebooks/tabular_examples/tree_based_models/Force%20Plot%20Colors.html?highlight=force%20plot
# https://www.sciencedirect.com/science/article/pii/S1532046418301758

"""
# From config (TODO, might want these all in rule)
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

iterations = int(snakemake.params["iterations"])

df_input = snakemake.input["in_file"]
out_file = snakemake.output["out_file"]
dataset = snakemake.params["dataset"]

# set graphics dimensions
_width = 9
_height = 6

# TODO review this
join_flag = True

out_prefix = snakemake.config["out"] + snakemake.config["path_" + dataset] + dataset

plot_list = []
plot_file_name_list = []


# saves to log file; does not print to screen
class Logger(object):
    def __init__(self):
        self.log = open(out_prefix + "-log.txt", "a")

    def write(self, message):
        self.log.write(message)

    def flush(self):
        pass


# placeholder dataframe for creating figure
def no_features_selected():
    _df = pd.DataFrame({"important_features": ["no_selected_features"],
                        "decoded_features": ["NA"],
                        "feature_importance_vals": [-100]})
    _df.to_csv(out_prefix + "-boruta-important.txt", index=False, sep="\t")


def fail_analysis(out_file):
    with open(out_file, "w") as file:
        file.write("Failed Analysis")
    no_features_selected()
    return


def percent_var_statement(step, oob):
    if step == "full_forest":
        # if first percent var explained is less than 5% then terminate
        print("\n% Variation Explained:")
        if float(oob) < 5.0:
            print("Out-of-bag (OOB) score: " + oob + "% \n\n*** Analysis " +
                  "failed due to low percent variation explained")
            fail_analysis(out_file)
            return

        print("Out-of-bag (OOB) score: " + oob + "%")

    if step == "shap_rerun":
        print("Out-of-bag (OOB) score (% variation explained) excluding "
              + "Boruta rejected features only for" +
              "visualization/interpretation): " + oob + "%")


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

    # For NSHAP Study
    # df = df.astype(str)
    # numeric_cols = ["weight_sel", "drugs_count", "np_count",
    #                 "weight_adj", "stratum", "cluster", "gender", 
    #                 "age", "happy"]
    # df[numeric_cols] = df[numeric_cols].astype(float)
    # END NSHAP
    numeric_column_list = df._get_numeric_data().columns.values
    categoric_columns = [i for i in list(df.columns) if
                         i not in numeric_column_list]

    if include_time == "yes":
        do_not_encode = [random_effect, sample_id, response_var]
    if include_time == "no":
        do_not_encode = [random_effect, sample_id, response_var,
                         "timepoint_explana"]

    study_n = str(df[random_effect].nunique())

    print("\nN (" + random_effect + "): ", study_n)

    # encode variables to improve looking at specific values
    vars_to_encode = [i for i in categoric_columns if i not in do_not_encode]
    df = dummy_var_creation(vars_to_encode, df, join_flag)

    # if response var is categoric convert to integers
    if df[response_var].dtype == "object":

        df = df.sort_values(by=response_var)

        # remove vars before encoding, random effect, response, sample ID
        df[response_var], u_maps = pd.factorize(df[response_var])

        print("\nNumeric mapping was created for categoric response" +
              "variable," + str(response_var) +
              ", (sorted alphabetically and factorized):")
        for i, category in enumerate(u_maps):
            print(f"-> {category} was converted to {i}")

        print("\n*** Percent variation explained is not optimal for " +
              "categorical response variables. \n*** Use caution with " +
              "interpretation.\n")

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

    x.to_csv(out_prefix + "-input-model-df.txt", index=False, sep="\t")
    df_features_in = pd.DataFrame(x.columns, columns=["input_features"])

    return x, z, c, y, df_features_in, study_n


def train_rf_model(x, y, rf_regressor, step):
    forest = rf_regressor.fit(x, y)
    oob = str(round(forest.oob_score_*100, 1))  # percent variation
    percent_var_statement(step, oob)
    plt.close("all")
    return forest, oob


def graphic_handling(plt, plot_file_name, p_title, width, 
                     height, tick_size, title_size):
    plot_file_name_list.append(plot_file_name)
    plt.title(p_title, fontdict={"fontsize": title_size})
    plt.tick_params(axis="both", labelsize=tick_size)
    p = plt.gcf()
    p.set_size_inches(width, height)
    plt.close("all")
    plot_list.append(p)


def list_split(input_list, group_size):
    max_items = group_size
    return [input_list[i:i+max_items]
            for i in range(0, len(input_list), max_items)]


def my_shap_summary_plots(feature_index_list, shap_values, x, list_num,
                          group_size):
    t = response_var + " values explained with selected features"
    plt_s = shap_values[:, feature_index_list]
    plt_x = x.iloc[:, feature_index_list]
    shap.summary_plot(shap_values=plt_s, features=plt_x, show=False,
                      sort=True, max_display=group_size)
    graphic_handling(plt, "-accepted-SHAP-summary-beeswarm-" +
                     str(list_num), t, _width, _height, 7, 10)
    shap.summary_plot(shap_values=plt_s, features=plt_x, show=False, sort=True,
                      plot_type="bar", max_display=group_size)
    graphic_handling(plt, "-accepted-SHAP-summary-bar-" + str(list_num), t,
                     _width, _height, 7, 10)


def shap_plots(shap_values, x, boruta_accepted, feature_importance):

    # get indices from Boruta for accepted features
    index_list_accepted = list(feature_importance[feature_importance[
        "col_name"].isin(boruta_accepted)].index)

    # Number of important features can vary, so split up graphs so
    # graphs aren't cluttered (note scales will vary)
    group_size = 10
    lists_to_graph = list_split(index_list_accepted, group_size=group_size)

    shap_values_df = pd.DataFrame(shap_values, columns=list(x.columns.values))
    shap_values_df.to_csv(out_prefix + "-SHAP-values-df.txt",
                          index=False, sep="\t")
    # graph each list of important features after
    for i in range(0, len(lists_to_graph)):
        my_shap_summary_plots(lists_to_graph[i], shap_values, x, i, group_size)

    col_count = 1
    # TODO save the dataframe and get vals 1
    for col in boruta_accepted:
        shap.dependence_plot(col, shap_values, x, show=False, x_jitter=0.03,
                             interaction_index=None)
        graphic_handling(plt, "-SHAP-dependence-plot-interaction-" +
                         str(col_count) + "-" + str(col),
                         "Dependence plot showing feature impact on " +
                         response_var, _width, _height, 7, 10)
        col_count += 1


def run_shap(trained_forest_model, x):
    shap_explainer = shap.TreeExplainer(trained_forest_model)
    shap_values = shap_explainer.shap_values(x)

    feature_names = list(x.columns.values)
    vals = np.abs(shap_values).mean(0)
    shap_imp = \
        pd.DataFrame(list(zip(feature_names, vals)), columns=["col_name",
                     "feature_importance_vals"])
    shap_imp.sort_values(by=["feature_importance_vals"],
                         ascending=False, inplace=True)
    shap_imp["feature_importance_vals"] = \
        shap_imp["feature_importance_vals"].apply(lambda x: round(x, 5))

    return shap_imp, shap_values


# TODO check percentile
def run_boruta_shap(forest, x, y):
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
                         boruta_width, _height, 6, 10)
    return feature_selector


def train_merf_model(x, z, c, y, model, step):
    mrf = MERF(max_iterations=iterations, fixed_effects_model=model)
    mrf.fit(x, z, c, y)
    forest = mrf.trained_fe_model
    oob = str(round(forest.oob_score_*100, 1))  # percent variation
    percent_var_statement(step, oob)
    return forest, oob


def main(df_input, out_file, random_effect, sample_id, response_var,
         join_flag):
    # make a log for printed information and a dataframe
    sys.stdout = Logger()

    # for log
    model_evaluation = "NA"
    var_explained = "NA"
    study_n = "NA"
    sample_n = "NA"
    input_features = "NA"
    accepted_features = "NA"
    tentative_features = "NA"
    rejected_features = "NA"

    log_df = {
        "Data": ["Model Evaluation", "% Variance Explained", "N (Study ID)",
                 "N (Samples)", "Input features", "Accepted Features",
                 "Tentative Features", "Rejected Features"],
        str.capitalize(dataset): ["NA", "NA", "NA", "NA", "NA", "NA", "NA",
                                  "NA"]
    }
    log_df = pd.DataFrame(log_df)

    df = pd.read_csv(df_input, sep="\t", na_filter=False)

    # if the input dataframe is empty; rule terminates
    if "Failed Analysis" in df.columns:
        model_evaluation = "Not performed"
        fail_analysis(out_file)

    else:
        study_n = str(df[random_effect].nunique())
        sample_n = str(df[sample_id].nunique())

        if len(np.unique(df[random_effect])) != len(df[random_effect]):
            fe_model_needed = False
            print("MODEL TYPE: Mixed Effects Random Forest (MERF) Regressor ")
        else:
            print("MODEL TYPE: Random Forest (RF) Regressor")
            fe_model_needed = True

        x, z, c, y, df_features_in, study_n = \
            setup_df_do_encoding(df=df, random_effect=random_effect,
                                 response_var=response_var,
                                 join_flag=join_flag, sample_id=sample_id)

        # Set up RF regressor for both fixed or mixed models
        def instantiate_rf():
            rf_regressor = \
                RandomForestRegressor(n_jobs=-1, n_estimators=n_estimators,
                                      oob_score=True,
                                      max_features=max_features,
                                      max_depth=7)  # for BorutaSHAP complexity
            return (rf_regressor)

        # Run either mixed effects or fixed effects RF models
        if not fe_model_needed:
            forest, oob = train_merf_model(x=x, z=z, c=c, y=y,
                                           model=instantiate_rf(),
                                           step="full_forest")
            plt.close("all")

        elif fe_model_needed:
            forest, oob = train_rf_model(x, y, instantiate_rf(),
                                         step="full_forest")

        # if first time running RF is < 5.0 percent variation, end analysis
        var_explained = str(oob) + "%"
        if float(oob) < 5.0:
            model_evaluation = "Fail: low variance"
            fail_analysis(out_file)

        else:
            # BorutaSHAP to select features
            boruta_explainer = run_boruta_shap(forest, x, y)

            # refit forest after excluding rejected features prior to
            # visualizing data. Essentially if data is visualized in SHAP plots
            # with all the rejected features, they're very hard to interpret
            # (flooded with meaningless info)
            # z, c, and y are the same as before

            x = x.drop(boruta_explainer.rejected, axis=1)

            if not fe_model_needed:
                forest, _ = train_merf_model(x=x, z=z, c=c, y=y,
                                             model=instantiate_rf(),
                                             step="shap_rerun")
            elif fe_model_needed:
                forest, _ = train_rf_model(x, y, instantiate_rf(),
                                           step="shap_rerun")

            shap_imp, shap_values = run_shap(forest, x)
            shap_plots(shap_values, x, boruta_explainer.accepted, shap_imp)

            shap_imp.to_csv(out_prefix + "SHAP-feature-imp.txt", sep='\t',
                            mode='a')

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

            df_boruta.to_csv(out_prefix + "-boruta-all-features.txt",
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
            df_features_in.to_csv(out_prefix + "-input-features.txt", sep="\t",
                                  index=False)

            if df_boruta.empty:
                no_features_selected()
            else:
                df_boruta.to_csv(out_prefix + "-boruta-important.txt",
                                 index=False, sep="\t")

            # logging data
            input_features = str(len(boruta_explainer.all_columns))
            accepted_features = str(len(boruta_explainer.accepted))
            tentative_features = str(len(boruta_explainer.tentative))
            rejected_features = str(len(boruta_explainer.rejected))
            print("\nBorutaSHAP features")
            print("Total input features: " + input_features)
            print("Accepted: " + accepted_features)
            print("Tentative: " + tentative_features)
            print("Rejected: " + rejected_features)
            plt.close("all")
            _build_result_pdf(out_file, out_prefix, plot_list,
                              plot_file_name_list)
            model_evaluation = "Pass"

    log_df[str.capitalize(dataset)] = [model_evaluation, var_explained,
                                       study_n, sample_n, input_features,
                                       accepted_features, tentative_features,
                                       rejected_features]
    log_df.to_csv(out_prefix + "-log-df.txt", index=False, sep="\t")

    out_path = snakemake.config["out"] + snakemake.config["path_" + dataset]
    files_in_folder = os.listdir(out_path)

    # remove some svgs because they're saved in PDF
    for file_name in files_in_folder:
        for n in ['SHAP-dependence-plot', 'features.svg']:
            if n in file_name:
                fp = os.path.join(out_path, file_name)
                os.remove(fp)


main(df_input=df_input, out_file=out_file, random_effect=random_effect,
     sample_id=sample_id, response_var=response_var, join_flag=join_flag)
