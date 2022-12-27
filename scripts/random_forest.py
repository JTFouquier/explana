

# TODO organize imports; these are a mess
import os

import shap
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
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

df_input = snakemake.input["in_file"]
pdf_report = snakemake.output["out_file"]

fe_or_me = snakemake.params["random_forest_type"]
random_effect = snakemake.params["random_effect"]
sample_id = snakemake.config["sample_id"]
response_var = snakemake.config["response_var"]
delta_flag = snakemake.params["delta_flag"]
iterations = int(snakemake.config["iterations"])

# set graphics dimensions
_width = 14
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


def setup_df_do_encoding(df, random_effect, response_var,
                         delta_flag, join_flag, sample_id):

    numeric_column_list = df._get_numeric_data().columns.values
    categoric_columns = [i for i in list(df.columns) if
                         i not in numeric_column_list]
    # exclude random_effect and sample_id
    vars_to_encode = \
        [i for i in categoric_columns if i not in [random_effect, sample_id]]

    df_dummy, dummy_dict = _dummy_var_creation(vars_to_encode, df, join_flag)
    # remove vars before encoding, random effect, response, sample ID
    df_encoded = df_dummy.drop([sample_id] + [random_effect] + vars_to_encode
                               + [response_var], axis=1)
    feature_list = [x for x in df_encoded.columns if "_reference" not in x]
    x = df_dummy[feature_list]  # Raw analysis and ALL (new, remove ref)
    df_final_input_features = pd.DataFrame(x.columns,
                                           columns=['InputFeatures'])
    df_final_input_features.to_csv(out_file_prefix + "-input-features.txt",
                                   index=False, sep="\t")

    z = np.ones((len(x), 1))
    c = df_dummy[random_effect]
    y = df_dummy[response_var]
    return x, z, c, y, feature_list, dummy_dict


def train_rf_model(x, y, rf_regressor):
    rf_regressor.fit(x, y)
    plt.plot([1, 2, 3, 4])
    plt.ylabel('place holder plot')
    training_stats_plot = plt.gcf()
    plt.close()
    return rf_regressor, training_stats_plot


def graphic_handling(plt, plot_file_name, p_title, width=_width, height=_height):
    plot_file_name_list.append(plot_file_name)
    plt.title(p_title, fontdict={"fontsize": 18})
    plt.tick_params(axis="both", labelsize=7)
    p = plt.gcf()
    p.set_size_inches(width, height)
    plt.close()
    plot_list.append(p)


def shap_plots(shap_values, x, boruta_accepted, boruta_accepted_tentative,
               feature_importance):

    # get indexes from Boruta for accepted and accepted+tentative features
    index_list_accepted = list(feature_importance[feature_importance[
        "col_name"].isin(boruta_accepted)].index)
    index_list_accepted_tentative = list(feature_importance[feature_importance[
        "col_name"].isin(boruta_accepted_tentative)].index)

    def my_shap_summary_plots(feature_index_list, shap_values, x,
                              accepted_flag):
        if accepted_flag:
            title_name_bee = "Accepted important features displayed with " \
                             "SHAP beeswarm plot"
            title_name_bar = "Accepted important features displayed with " \
                             "SHAP bar plot"
            file_name_bee = "-accepted-SHAP-summary-beeswarm"
            file_name_bar = "-accepted-SHAP-summary-bar"
        if not accepted_flag:
            title_name_bee = "Accepted and tentative important features " \
                             "displayed with SHAP beeswarm plot"
            title_name_bar = "Accepted and tentative important features " \
                             "displayed with SHAP bar plot"
            file_name_bee = "-accepted-and-tentative-SHAP-summary-beeswarm"
            file_name_bar = "-accepted-and-tentative-SHAP-summary-bar"
        shap.summary_plot(shap_values=shap_values[:, feature_index_list],
                          features=x.iloc[:, feature_index_list], show=False,
                          plot_size=(_width, _height), sort=False)
        graphic_handling(plt, file_name_bee, title_name_bee)
        shap.summary_plot(shap_values=shap_values[:, feature_index_list],
                          features=x.iloc[:, feature_index_list], show=False,
                          plot_size=(_width, _height), sort=False,
                          plot_type="bar")
        graphic_handling(plt, file_name_bar, title_name_bar)

    my_shap_summary_plots(index_list_accepted, shap_values, x, True)
    # Uncomment to print accepted and tentative features - not necessary
    # my_shap_summary_plots(index_list_accepted_tentative, shap_values,
    # x, False)

    col_count = 1
    # TODO make either accepted or accepted and tentative
    for col in boruta_accepted:
        shap.dependence_plot(col, shap_values, x, show=False, x_jitter=0.03)
        p_title = "Dependence plot for accepted features from BorutaShap \n" \
                  "displayed with SHAP values (with interaction)"
        graphic_handling(plt, "-dependence-plot-interaction-" +
                         str(col_count) + str(col), p_title)


def run_shap(trained_forest_model, x):
    explainer = shap.TreeExplainer(trained_forest_model)
    shap_values = explainer.shap_values(x)

    feature_names = list(x.columns.values)
    vals = np.abs(shap_values).mean(0)
    shap_imp = \
        pd.DataFrame(list(zip(feature_names, vals)), columns=['col_name',
                     'feature_importance_vals'])
    shap_imp.sort_values(by=['feature_importance_vals'],
                         ascending=False, inplace=True)
    shap_imp['feature_importance_vals'] = \
        shap_imp['feature_importance_vals'].apply(lambda x: round(x, 5))

    return shap_imp, shap_values


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
        # print(feature_group)
        # print(boruta_dict)
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

def main(df_input, out_file, random_forest_type, random_effect, sample_id,
         response_var, delta_flag, join_flag):
    df = pd.read_csv(df_input, sep="\t", na_filter=False)

    x, z, c, y, feature_list, dummy_dict = \
        setup_df_do_encoding(df=df, random_effect=random_effect,
                             response_var=response_var, delta_flag=delta_flag,
                             join_flag=join_flag, sample_id=sample_id)

    rf_regressor = \
        RandomForestRegressor(n_jobs=-1, n_estimators=n_estimators,
                              oob_score=True, max_features=0.90)
    rf_param_dict = rf_regressor.get_params()
    print("RF ESTIMATORS (TREES): ", rf_param_dict['n_estimators'])
    # print("Random Forest Regressor parameters:\n" +
    #       str(rf_regressor.get_params()))


    # TODO!!!!  if two timepoints original MERF deltas RF
    if random_forest_type == "mixed":
        forest, mrf = \
            run_mixed_effects_random_forest(x=x, z=z, c=c, y=y,
                                            model=rf_regressor,
                                            step="full_forest")

        plot_merf_training_stats(mrf, 12)  # TODO fix
        plt.tight_layout()
        training_stats_plot = plt.gcf()
        plt.close()

    else:
        # TODO
        mrf, training_stats_plot = train_rf_model(x, y)

    # BorutaSHAP to select features
    feature_selector = run_boruta_shap(forest, x, y)
    # print("feature shap values")
    # print(feature_selector.shap_values)
    # print("feature accepted")
    # print(feature_selector.accepted)
    # print("feature accepted columns")
    # print(feature_selector.accepted_columns)
    # SHAP to explain/detail features' values' effects on response variable
    # With rerunning SHAP on only important features (i.e. adding all
    # unimportant features causes unrealistic plots and poor interpretation)

    # refit forest after excluding rejected features prior to visualizing data
    # Essentially if data is visualized in SHAP plots with all the rejected
    # features, they're very hard to interpret (flooded with meaningless info)
    refit = "True"
    if refit == "True":
        # z, c, and y are the same as before
        x = x.drop(feature_selector.rejected, axis=1)
        forest, mrf = \
            run_mixed_effects_random_forest(x=x, z=z, c=c, y=y,
                                            model=rf_regressor,
                                            step="shap_rerun")
    else:
        print("SHAP plots visualized using all input features")

    shap_imp, shap_values = run_shap(forest, x)
    shap_plots(shap_values, x, feature_selector.accepted,
               feature_selector.accepted + feature_selector.tentative,
               shap_imp)

    shap_imp.to_csv(out_file_prefix + "SHAP-feature-imp.txt",
                    sep='\t', mode='a')

    # find prefix for encoded values "ENC_Location_is_1_3" decoded is Location
    # TODO find key words that can't be in column names (_is_, ENC_, etc)
    decoded_top_features_list = []
    for i in feature_selector.accepted:
        # print("feature in accepted")
        # print(i)
        i_for_dict = i.split("_is_")[0] + "_is"  # TODO important word features
        try:
            decoded = dummy_dict[i_for_dict]
        except KeyError:
            decoded = i
        decoded_top_features_list.append(''.join(decoded))

    df_boruta = pd.DataFrame({'important.features': feature_selector.accepted,
                              'decoded.features': decoded_top_features_list})
    df_boruta.to_csv(out_file_prefix + "-boruta-important.txt", index=False,
                     sep="\t")
    print("\nBorutaSHAP features")
    print("Total input features: " +
          str(len(feature_selector.all_columns)))
    print("Accepted: " + str(len(feature_selector.accepted)))
    print("Tentative: " + str(len(feature_selector.tentative)))
    print("Accepted and Tentative: " + str(len(feature_selector.accepted) +
                                         len(feature_selector.tentative)))
    print("Rejected: " + str(len(feature_selector.rejected)))
    print("\nTop BorutaSHAP Accepted Features (max 30 displayed)\n")
    print(df_boruta.head(30))
    plot_partial_dependence(forest, x, features=feature_selector.accepted)
    plot_some_partial_dependence = plt.gcf()
    plot_some_partial_dependence.set_size_inches(_width * 1.4, _height * 1.4)
    plot_list.append(plot_some_partial_dependence)
    plot_file_name_list.append("-partial-dependence")
    plot_list.append(training_stats_plot)
    plot_file_name_list.append("-training-stats")
    _build_result_pdf(out_file, out_file_prefix, plot_list, plot_file_name_list)


main(df_input=df_input, out_file=pdf_report, random_forest_type=fe_or_me,
     random_effect=random_effect, sample_id=sample_id,
     response_var=response_var, delta_flag=delta_flag, join_flag=join_flag)
