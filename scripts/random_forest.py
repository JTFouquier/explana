

# TODO organize imports; these are a mess
import os

import shap
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import re
import math
import yaml
from sklearn.inspection import plot_partial_dependence

from merf.merf import MERF
from merf.viz import plot_merf_training_stats
from BorutaShap import BorutaShap
# from boruta_shap_mine import BorutaShap
from sklearn.ensemble import RandomForestRegressor

# TODO switch to OneHotEncoder
from hot_encoding import _dummy_var_creation
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
n_estimators = int(snakemake.config["n_estimators"])
borutaSHAP_trials = int(snakemake.config["borutaSHAP_trials"])
borutaSHAP_percentile = int(snakemake.config["borutaSHAP_percentile"])
enc_percent_threshold = float(snakemake.config["enc_percent_threshold"])

df_input = snakemake.input["in_file"]
pdf_report = snakemake.output["out_file"]

dataset = snakemake.params["dataset"]

random_effect = snakemake.config["random_effect"]
sample_id = snakemake.config["sample_id"]
response_var = snakemake.config["response_var"]
iterations = int(snakemake.config["iterations"])
# set graphics dimensions
_width = 13
_height = 8

# TODO review this
join_flag = True
out_file_prefix = pdf_report.split(".pdf")[0]

plot_list = []
plot_file_name_list = []

# TODO this is not ideal for logging, but it works (issues with stream/file)
class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open(out_file_prefix + "-log.txt", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        pass

sys.stdout = Logger()


def setup_df_do_encoding(df, random_effect, response_var, join_flag, sample_id,
                         fe_model_needed):

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

    # whether or not fixed effects or mixed effects model is needed
    if not fe_model_needed:
        drop_list = [random_effect, sample_id]
    if fe_model_needed:
        drop_list = [sample_id]
    print("Input features before encoding: ", df.shape[1])

    # encode variables to improve looking at specific values
    vars_to_encode = [i for i in categoric_columns if i not in drop_list]
    df_dummy = _dummy_var_creation(vars_to_encode, df, join_flag)
    print("Input features after encoding: ", df_dummy.shape[1])

    # remove vars before encoding, random effect, response, sample ID
    df_encoded = df_dummy.drop(drop_list + vars_to_encode + [response_var], axis=1)

    def remove_low_percentage_encoded_vars(df_encoded, percent):
        # TODO make function for removing columns
        # yes/1 counts per column; drop columns in 
        df_complete = df_encoded
        non_encoded_cols = [col for col in df_encoded.columns if not 'ENC_' in col]
        df_categoric = df_encoded.drop(non_encoded_cols, axis=1)

        value_counts = df_categoric.eq(1).sum()
        target_value = round(len(df_categoric)*(percent/100))

        # list of column names to drop from df
        columns_to_drop = list(value_counts[value_counts < target_value].index)
        df = df_complete.drop(columns_to_drop, axis=1)
        return(df)

    if enc_percent_threshold >= 0:
        df_encoded = \
        remove_low_percentage_encoded_vars(df_encoded=df_encoded, 
                                           percent=enc_percent_threshold)
        
    print("Input features after dropping vars below " + 
          str(enc_percent_threshold) + " percent: ", 
          df_encoded.shape[1])
    
    feature_list = [x for x in df_encoded.columns]

    x = df_dummy[feature_list]  # Raw analysis and ALL (new, remove ref)
    z = np.ones((len(x), 1))
    c = df_dummy[random_effect]
    y = df_dummy[response_var]

    x.to_csv(out_file_prefix + "-input-model-df.txt", 
             index=False, sep="\t")

    df_final_input_features = pd.DataFrame(x.columns,
                                           columns=['InputFeatures'])
    df_final_input_features.to_csv(out_file_prefix + "-input-features.txt",
                                   index=False, sep="\t")

    return x, z, c, y, feature_list

# TODO add rerun step here
def train_rf_model(x, y, rf_regressor):
    forest = rf_regressor.fit(x, y)
    plt.plot([1, 2, 3, 4])
    plt.ylabel('place holder plot')
    training_stats_plot = plt.gcf()
    # plt.close()
    plt.close('all')
    return forest, training_stats_plot


def graphic_handling(plt, plot_file_name, p_title, width=_width, height=_height):
    plot_file_name_list.append(plot_file_name)
    plt.title(p_title, fontdict={"fontsize": 10})
    plt.tick_params(axis="both", labelsize=7)
    p = plt.gcf()
    p.set_size_inches(width, height)
    plt.close('all')
    plot_list.append(p)


def list_split(input_list, group_size):
    max_items = group_size
    return [input_list[i:i+max_items]
            for i in range(0, len(input_list), max_items)]

def shap_plots(shap_values, x, boruta_accepted, feature_importance, 
               shap_explainer):

    # get indices from Boruta for accepted features
    index_list_accepted = list(feature_importance[feature_importance[
        "col_name"].isin(boruta_accepted)].index)

    # Number of important features can vary, so split up graphs so 
    # graphs aren't cluttered (note scales will vary)
    group_size = 10
    lists_to_graph = list_split(index_list_accepted, group_size = group_size)

    def my_shap_summary_plots(feature_index_list, shap_values, x,
                              accepted_flag, list_num):
        if accepted_flag:
            title_name_bee = "Accepted important features displayed with " \
                             "SHAP beeswarm plot"
            title_name_bar = "Accepted important features displayed with " \
                             "SHAP bar plot"
        plt_s = shap_values[:, feature_index_list]
        plt_x = x.iloc[:, feature_index_list]

        shap.summary_plot(shap_values=plt_s, features=plt_x, show=False, sort=True, 
                          max_display=group_size)
        graphic_handling(plt, "-accepted-SHAP-summary-beeswarm-" + str(list_num), 
                         title_name_bee)
        shap.summary_plot(shap_values=plt_s, features=plt_x, show=False, sort=True,
                          plot_type="bar", max_display=group_size)
        graphic_handling(plt, "-accepted-SHAP-summary-bar-" + str(list_num), 
                         title_name_bar)

    # save SHAP values to inspect TODO remove later
    shap_values_df = pd.DataFrame(shap_values, columns = list(x.columns.values))
    shap_values_df.to_csv(out_file_prefix + "-SHAP-values-df.txt",
                          index=False, sep="\t")
    # graph each list of important features after 
    for l in range(0, len(lists_to_graph)):
        my_shap_summary_plots(lists_to_graph[l], shap_values, x, True, l)

    col_count = 1
    ### TODO save the dataframe and get vals 1
    for col in boruta_accepted:
        shap.dependence_plot(col, shap_values, x, show=False, x_jitter=0.03)
        graphic_handling(plt, "-SHAP-dependence-plot-interaction-" +
                         str(col_count) + "-" + str(col), 
                         p_title="Dependence plot for accepted features " +
                         "from BorutaSHAP \n displayed with SHAP values")
        col_count += 1


def run_shap(trained_forest_model, x):
    shap_explainer = shap.TreeExplainer(trained_forest_model)
    shap_values = shap_explainer.shap_values(x)

    feature_names = list(x.columns.values)
    vals = np.abs(shap_values).mean(0)
    shap_imp = \
        pd.DataFrame(list(zip(feature_names, vals)), columns=['col_name',
                     'feature_importance_vals'])
    shap_imp.sort_values(by=['feature_importance_vals'],
                         ascending=False, inplace=True)
    shap_imp['feature_importance_vals'] = \
        shap_imp['feature_importance_vals'].apply(lambda x: round(x, 5))

    return shap_imp, shap_values, shap_explainer


# TODO check percentile
def run_boruta_shap(forest, x, y):
    feature_selector = \
        BorutaShap(model=forest, importance_measure='shap',
                   classification=False, percentile=borutaSHAP_percentile,
                   pvalue=0.05)

    feature_selector.fit(X=x, y=y, n_trials=borutaSHAP_trials,
                         train_or_test="train", sample=False, verbose=False)
    for feature_group in ['all', 'accepted', 'tentative', 'rejected']:
        boruta_dict = \
            _boruta_shap_plot(feature_selector, which_features=feature_group,
                              figsize=(_width, _height), X_size=8)
        graphic_handling(plt, "-boruta-" + feature_group + "-features",
                         "Boruta " + feature_group + " features",
                         width=_width, height=_height)
    return feature_selector


def run_mixed_effects_random_forest(x, z, c, y, model, step):
    mrf = MERF(max_iterations=iterations, fixed_effects_model=model)
    mrf.fit(x, z, c, y)
    forest = mrf.trained_fe_model
    oob = str(round(forest.oob_score_*100, 1))  # percent variation
    if step == "full_forest":
        print("\nOut-of-bag (OOB) score "
              "(% variation explained): " + oob + "%")
    if step == "shap_rerun":
        print("Out-of-bag (OOB) score "
              "(% variation explained excluding Boruta rejected "
              "features only for visualization/interpretation): " + oob + "%")
    return forest, mrf


def main(df_input, out_file, random_effect, sample_id, response_var,
         join_flag):
    df = pd.read_csv(df_input, sep="\t", na_filter=False)

    if len(np.unique(df[random_effect])) != len(df[random_effect]):
        fe_model_needed = False
        print("MODEL TYPE: Mixed Effects Random Forest (MERF) Regressor ")
    else:
        print("MODEL TYPE: Random Forest (RF) Regressor")
        fe_model_needed = True

    x, z, c, y, feature_list = \
        setup_df_do_encoding(df=df, random_effect=random_effect,
                             response_var=response_var,
                             join_flag=join_flag, sample_id=sample_id,
                             fe_model_needed=fe_model_needed)

    # Set up RF regressor for both fixed or mixed models
    def instantiate_rf(num_features):
        rf_regressor = \
            RandomForestRegressor(n_jobs=-1, 
                                  n_estimators=n_estimators,
                                  oob_score=True,
                                  max_depth=round(math.sqrt(num_features)))
        return(rf_regressor)

    # Run either mixed effects or fixed effects RF models
    if not fe_model_needed:
        forest, mrf = \
            run_mixed_effects_random_forest(x=x, z=z, c=c, y=y,
                                            model=instantiate_rf(len(x.columns)),
                                            step="full_forest")
        plot_merf_training_stats(mrf, 12)  # TODO fix
        plt.tight_layout()
        training_stats_plot = plt.gcf()
        plt.close('all')

    elif fe_model_needed:
        forest, training_stats_plot = \
        train_rf_model(x, y, instantiate_rf(len(x.columns)))

    # BorutaSHAP to select features
    boruta_explainer = run_boruta_shap(forest, x, y)

    # refit forest after excluding rejected features prior to visualizing data
    # Essentially if data is visualized in SHAP plots with all the rejected
    # features, they're very hard to interpret (flooded with meaningless info)
    # z, c, and y are the same as before
    x = x.drop(boruta_explainer.rejected, axis=1)

    if not fe_model_needed:
        forest, mrf = \
        run_mixed_effects_random_forest(x=x, z=z, c=c, y=y,
                                        model=instantiate_rf(len(x.columns)),
                                        step="shap_rerun")
    elif fe_model_needed:
        forest, _ = train_rf_model(x, y, instantiate_rf(len(x.columns)))

    shap_imp, shap_values, shap_explainer = run_shap(forest, x)
    shap_plots(shap_values, x, boruta_explainer.accepted, shap_imp, 
               shap_explainer)

    shap_imp.to_csv(out_file_prefix + "SHAP-feature-imp.txt",
                    sep='\t', mode='a')

    # find prefix for encoded values "ENC_Location_is_1_3" decoded is Location
    # TODO find key words that can't be in column names (_is_, ENC_, etc)
    decoded_top_features_list = []
    shap_imp_score_list = []
    # get original column name if feature was encoded; else use feature name
    for feature_name in boruta_explainer.accepted:
        try:
            r = re.search('ENC_(.*)_is_', feature_name)
            decoded = r.group(1)
        except AttributeError:
            decoded = feature_name

        decoded_top_features_list.append(''.join(decoded))
        shap_imp_score = shap_imp.loc[shap_imp['col_name'] == feature_name,
                                      'feature_importance_vals']
        shap_imp_score_list.append(shap_imp_score.iloc[0])

    df_boruta = pd.DataFrame({'important.features': boruta_explainer.accepted,
                              'decoded.features': decoded_top_features_list,
                              'feature_importance_vals': shap_imp_score_list})
    df_boruta = df_boruta.sort_values(by=['feature_importance_vals'],
                                      ascending=False)
    df_boruta.to_csv(out_file_prefix + "-boruta-important.txt", index=False,
                     sep="\t")
    print("\nBorutaSHAP features")
    print("Total input features: " + str(len(boruta_explainer.all_columns)))
    print("Accepted: " + str(len(boruta_explainer.accepted)))
    print("Tentative: " + str(len(boruta_explainer.tentative)))
    print("Accepted and Tentative: " + str(len(boruta_explainer.accepted) +
                                           len(boruta_explainer.tentative)))
    print("Rejected: " + str(len(boruta_explainer.rejected)))
    plot_partial_dependence(forest, x, features=boruta_explainer.accepted)
    plot_some_partial_dependence = plt.gcf()
    plot_some_partial_dependence.set_size_inches(_width * 1.4, _height * 1.4)
    plot_list.append(plot_some_partial_dependence)
    plot_file_name_list.append("-partial-dependence")
    plot_list.append(training_stats_plot)
    plot_file_name_list.append("-training-stats")
    plt.close('all')
    _build_result_pdf(out_file, out_file_prefix, plot_list, plot_file_name_list)


main(df_input=df_input, out_file=pdf_report, random_effect=random_effect,
     sample_id=sample_id, response_var=response_var, join_flag=join_flag)

