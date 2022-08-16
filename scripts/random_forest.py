

# TODO organize imports; these are a mess
import shap
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from sklearn.inspection import plot_partial_dependence

from merf.merf import MERF
from merf.viz import plot_merf_training_stats
from BorutaShap import BorutaShap

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

df_deltas = snakemake.input["in_file"]
pdf_report = snakemake.output["out_file"]

fe_or_me = snakemake.params["random_forest_type"]
random_effect = snakemake.params["random_effect"]
sample_id = snakemake.config["sample_id"]
response_var = snakemake.config["response_var"]
delta_flag = snakemake.params["delta_flag"]
iterations = int(snakemake.config["iterations"])
re_timepoint = snakemake.params["re_timepoint"]

# set graphics dimensions
_width = 14
_height = 8
_font_title = 14

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


def setup_df_do_encoding(df, random_effect, vars_to_encode, response_var,
                         delta_flag, join_flag):
    df.to_csv(out_file_prefix + "-df-before-encoding.txt", index=False,
              sep="\t")

    df, dummy_dict = _dummy_var_creation(vars_to_encode, df, join_flag)
    df.to_csv(out_file_prefix + "-df-after-encoding.txt", index=False,
              sep="\t")

    drop = [sample_id] + \
           [random_effect] + vars_to_encode + [response_var]

    df_encoded_cols_dropped = df.drop(drop, axis=1)

    feature_list = df_encoded_cols_dropped.columns
    feature_list = [x for x in feature_list if "_reference" not in x]
    feature_list = [x for x in feature_list if "StudyID" not in x]
    print("Input features after encoding:\n")
    print(feature_list)
    x = df[feature_list]  # Raw analysis and ALL (new, remove ref)
    x.to_csv(out_file_prefix + "-input-features-to-model.txt", index=False,
             sep="\t")

    if re_timepoint == "re_timepoint":  # TODO fix this
        z = np.ones((len(x), 1))
    if re_timepoint == "no_re":
        z = np.ones((len(x), 1))
    clusters = df["StudyID"]
    y = df[response_var]
    x.to_csv(out_file_prefix + "-x.txt", index=False, sep="\t")
    return x, z, clusters, y, feature_list, dummy_dict


def train_rf_model(x, y, rf_regressor):
    rf_regressor.fit(x, y)
    plt.plot([1, 2, 3, 4])
    plt.ylabel('place holder plot')
    training_stats_plot = plt.gcf()
    plt.close()
    return rf_regressor, training_stats_plot


def graphic_handling(plt, plot_file_name, width=_width, height=_height):
    plot_file_name_list.append(plot_file_name)
    plt.tight_layout()
    p = plt.gcf()
    p.set_size_inches(width, height)
    plt.close()
    plot_list.append(p)


def shap_plots(shap_values, x, boruta_accepted, boruta_accepted_tentative,
               feature_importance):

    index_list = list(feature_importance[feature_importance[
        "col_name"].isin(boruta_accepted)].index)
    # Do not sort because it's using BorutaShap rankings
    shap.summary_plot(shap_values=shap_values[:, index_list],
                      features=x.iloc[:, index_list], show=False,
                      plot_size=(_width, _height), sort=False)
    plt.title("Accepted important features from BorutaShap displayed "
              "with SHAP Summary Plot", fontdict={'fontsize': _font_title})
    graphic_handling(plt, "-accepted_SHAP")
    index_list = list(feature_importance[feature_importance[
        "col_name"].isin(boruta_accepted_tentative)].index)
    shap.summary_plot(shap_values=shap_values[:, index_list],
                      features=x.iloc[:, index_list], show=False,
                      plot_size=(_width, _height), sort=False)
    plt.title("Accepted and tentative important features from BorutaShap "
              "displayed using SHAP values",
              fontdict={'fontsize': _font_title})
    graphic_handling(plt, "-accepted-tentative-SHAP-summary")

    col_count = 1
    # TODO make either accepted or accepted and tentative
    for col in boruta_accepted:
        shap.dependence_plot(col, shap_values, x, show=False, x_jitter=0.03)
        plt.title("Dependence plot for accepted features from "
                  "BorutaShap \ndisplayed with SHAP values (with interaction)",
                  fontdict={'fontsize': _font_title})
        graphic_handling(plt, "-dependence-plot-interaction-" +
                         str(col_count) + str(col))


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
    boruta_dict = {'accepted': [], 'tentative': [],
                   'rejected': [], 'all': []}

    for feature_group in ['all', 'accepted', 'tentative', 'rejected']:
        boruta_dict[feature_group] = \
            _boruta_shap_plot(feature_selector, which_features=feature_group,
                              figsize=(_width, _height), X_size=12)
        graphic_handling(plt, "-boruta-" + feature_group + "-features",
                         width=_width, height=_height)
    return boruta_dict


def main(df_input, out_file, random_forest_type, random_effect, sample_id,
         response_var, delta_flag, join_flag):
    df = pd.read_csv(df_input, sep="\t")

    numeric_column_list = list(df._get_numeric_data().columns)
    column_list = list(df.columns)
    categoric_columns = [i for i in column_list if
                         i not in numeric_column_list]

    # exclude random effect cols and sample_id
    # TODO do I want this in model? Fixed vs random effects?
    to_drop = [random_effect, sample_id]
    encode_this_list = [i for i in categoric_columns if i not in to_drop]

    x, z, clusters, y, feature_list, dummy_dict = \
        setup_df_do_encoding(df=df, random_effect=random_effect,
                             vars_to_encode=encode_this_list,
                             response_var=response_var, delta_flag=delta_flag,
                             join_flag=join_flag)

    rf_regressor = RandomForestRegressor(n_jobs=-2, n_estimators=n_estimators,
                                         oob_score=True, max_features='auto')

    def run_mixed_effects_random_forest(x, z, clusters, y, model):
        mrf = MERF(max_iterations=iterations, fixed_effects_model=model)
        mrf.fit(x, z, clusters, y)
        forest = mrf.trained_fe_model

        return forest, mrf

    if random_forest_type == "mixed":

        forest, mrf = run_mixed_effects_random_forest(x, z, clusters, y,
                                                      rf_regressor)
        print("\nRandom Forest out-of-bag (OOB) score: " +
              str(round(forest.oob_score_, 3)))
        print("\n\nrf_regressor params:" + str(rf_regressor.get_params()))

        plot_merf_training_stats(mrf, 12)  # TODO fix
        plt.tight_layout()
        training_stats_plot = plt.gcf()
        plt.close()

    else:
        # TODO
        mrf, training_stats_plot = train_rf_model(x, y)

    # BorutaSHAP to select features
    boruta_dict = run_boruta_shap(forest, x, y)
    boruta_accepted = boruta_dict['accepted']
    boruta_accepted_tentative = \
        boruta_dict['accepted'] + boruta_dict['tentative']

    # SHAP to explain/detail features' values' effects on response variable
    shap_imp, shap_values = run_shap(forest, x)
    shap_plots(shap_values, x, boruta_accepted, boruta_accepted_tentative,
               shap_imp)

    shap_imp.to_csv(out_file_prefix + "SHAP-feature-imp.txt",
                    sep='\t', mode='a')

    # find prefix for encoded values "ENC_Location_is_1_3" decoded is Location
    # TODO find key words that can't be in column names (_is_, ENC_, etc)
    decoded_top_features_list = []
    for i in boruta_accepted:
        i_for_dict = i.split("_is_")[0] + "_is"  # TODO important word features
        try:
            decoded = dummy_dict[i_for_dict]
        except KeyError:
            decoded = i
        decoded_top_features_list.append(''.join(decoded))

    df_boruta = pd.DataFrame({'important.features': boruta_accepted,
                              'decoded.features': decoded_top_features_list})
    df_boruta.to_csv(out_file_prefix + "-boruta-important.txt", index=False,
                     sep="\t")
    print("\nRanking of top BorutaSHAP Features (max 20 displayed)\n")
    print(df_boruta.head(20))
    plot_partial_dependence(forest, x, features=boruta_accepted)
    plot_some_partial_dependence = plt.gcf()
    plot_some_partial_dependence.set_size_inches(_width * 1.4, _height * 1.4)
    plot_list.append(plot_some_partial_dependence)
    plot_file_name_list.append("-partial-dependence")
    plot_list.append(training_stats_plot)
    plot_file_name_list.append("-training-stats")
    _build_result_pdf(out_file, out_file_prefix, plot_list, plot_file_name_list)


main(df_input=df_deltas, out_file=pdf_report, random_forest_type=fe_or_me,
     random_effect=random_effect, sample_id=sample_id,
     response_var=response_var, delta_flag=delta_flag, join_flag=join_flag)
