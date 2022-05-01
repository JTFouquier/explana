

# TODO organize imports; these are a mess
import shap
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import logging
import sys
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.inspection import plot_partial_dependence

from merf.merf import MERF
from merf.viz import plot_merf_training_stats
from BorutaShap import BorutaShap

from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error

from sklearn.model_selection import train_test_split

# TODO switch to OneHotEncoder
from hot_encoding import _dummy_var_creation

from scipy.stats.stats import pearsonr

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

# access values from snakemake's Snakefile
in_file = snakemake.input["in_file"]
out_file = snakemake.output["out_file"]

random_forest_type = snakemake.params["random_forest_type"]
random_effect = snakemake.params["random_effect"]
sample_ID = snakemake.params["sample_ID"]
drop_rows = snakemake.params["drop_rows"]
subset_rows = snakemake.params["subset_rows"]
drop_cols = snakemake.params["drop_cols"]
subset_cols = snakemake.params["subset_cols"]
response_var = snakemake.params["response_var"]
delta_flag = snakemake.params["delta_flag"]
max_iters = snakemake.params["iterations"]
re_timepoint = snakemake.params["re_timepoint"]

# set graphics dimensions
width = 18
height = 10

# features to show
max_display = 8
# for partial dep plots and more
top_features = 8

# TODO review this
join_flag = True
out_file_prefix = out_file.split(".pdf")[0]


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


def setup_df_do_encoding(df, random_effect, vars_to_encode,
                         response_var, delta_flag, join_flag):
    df.to_csv(out_file_prefix + "-df-before-encoding.txt", index=False,
              sep="\t")

    df, dummy_dict = _dummy_var_creation(vars_to_encode, df,
                                         join_flag)
    df.to_csv(out_file_prefix + "-df-after-encoding.txt", index=False,
              sep="\t")

    if delta_flag == "deltas":
        drop = ["StudyID.Timepoint", "SampleID"] + \
               [random_effect] + vars_to_encode + [response_var]
    if delta_flag == "raw":
        drop = [random_effect] + vars_to_encode + [response_var]

    df_encoded_cols_dropped = df.drop(drop, axis=1)
    # TODO remove columns that have "ref"
    feature_list = df_encoded_cols_dropped.columns

    df_encoded_cols_dropped.to_csv(out_file_prefix +
                                   "-dummy-variables-dropped-columns.txt",
                                   index=False, sep="\t")

    feature_list = [x for x in feature_list if "_reference" not in x]

    x = df[feature_list]  # Raw analysis and ALL (new, remove ref)

    if re_timepoint == "re_timepoint": # TODO fix this
        z = np.ones((len(x), 1))
    if re_timepoint == "no_re":
        z = np.ones((len(x), 1))
    clusters = df["StudyID"]
    y = df[response_var]
    x.to_csv(out_file_prefix + "-x.txt", index=False, sep="\t")
    return x, z, clusters, y, feature_list, dummy_dict


def merf_training_stats_plot(mrf, num_clusters):
    # print training stats TODO understand this plot more
    plot_merf_training_stats(mrf, num_clusters_to_plot=num_clusters)
    training_stats_plot = plt.gcf()
    plt.close()
    return training_stats_plot


def scatter_plot(y_axis, x_axis, forest_score):
    y = y_axis
    x = x_axis
    x_axis = np.array(x_axis)
    y_axis = np.array(y_axis)
    plt.scatter(y_axis, x_axis)
    plt.annotate("r-squared = {:.3f}".format(forest_score),
                 (min(y), max(x)))
    corr_plot = plt.gcf()
    plt.close()
    return corr_plot


def train_merf_model(x, z, clusters, y, max_iters):
    # TODO CHANGE THIS later! NEEDS MORE ITERATIONS
    mrf = MERF(max_iterations=max_iters,
               fixed_effects_model=RandomForestRegressor(n_estimators=100,
                                                         n_jobs=-2,
                                                         oob_score=True,
                                                         max_features="auto"),)
    mrf.fit(x, z, clusters, y)
    training_stats_plot = plot_merf_training_stats(mrf, 12)
    plt.close()
    return mrf, training_stats_plot


def train_rf_model(x, y):
    mrf = RandomForestRegressor()
    mrf.fit(x, y)
    plt.plot([1, 2, 3, 4])
    plt.ylabel('place holder plot')
    training_stats_plot = plt.gcf()
    plt.close()

    return mrf, training_stats_plot


def graphic_handling(graphic):
    p = graphic
    plt.close()
    return p


def important_shapley_features(shap_values, x):
    feature_names = list(x.columns.values)
    vals = np.abs(shap_values).mean(0)
    feature_importance = pd.DataFrame(list(zip(feature_names, vals)),
                                      columns=['col_name',
                                               'feature_importance_vals'])
    feature_importance.sort_values(by=['feature_importance_vals'],
                                   ascending=False, inplace=True)
    feature_importance['feature_importance_vals'] = \
        feature_importance['feature_importance_vals'].apply(lambda x:
                                                            round(x, 4))
    return feature_importance


def shap_explainer(trained_forest_model, x):
    explainer = shap.TreeExplainer(trained_forest_model)
    shap_values = explainer.shap_values(x)
    plot_list = []

    # Bee swarm plot
    shap.summary_plot(shap_values=shap_values, features=x, show=False,
                      plot_size=(width, height), max_display=max_display)
    plot_list.append(graphic_handling(plt.gcf()))

    # Bar plot
    shap.summary_plot(shap_values=shap_values, features=x, plot_type="bar",
                      show=False, plot_size=(width, height),
                      max_display=max_display)
    plot_list.append(graphic_handling(plt.gcf()))

    # Feature Shapley values starting with best at explaining the model
    feature_importance = important_shapley_features(shap_values, x)

    top_features_list = list(feature_importance['col_name'][0:top_features])
    for col in top_features_list:
        shap.dependence_plot(col, shap_values, x, show=False)
        p = plt.gcf()
        p.set_size_inches(width, height)
        plot_list.append(p)
        plt.close()

    return plot_list, feature_importance, top_features_list


def build_result_pdf(out_file, plot_list):
    """
    Build PDF containing multiple .gcf() objects (get current figure)
    Args:
        pdf_name: PDF name for storing multiple plots
        plot_list: A list of .gcf() objects
    """
    pp = PdfPages(out_file)
    for i in plot_list:
        pp.savefig(i)
    pp.close()


def main(in_file, out_file, random_forest_type, random_effect, sample_ID,
         response_var, delta_flag, join_flag):
    df = pd.read_csv(in_file, sep="\t")

    # subsets of data
    for k, v in drop_rows.items():
        try:
            df = df.loc[(df[k] != v)]
        except:
            pass

    for k, v in subset_rows.items():
        try:
            df = df.loc[(df[k] == v)]
        except:
            pass

    # Can't have both
    if subset_cols:
        if drop_cols:
            print("Cannot have features (arguments) in 'drop_cols' and "
                  "'subset_cols' parameters")

    if subset_cols:
        df = df[subset_cols]
    else:
        pass

    df = df.drop(drop_cols, axis=1)

    numeric_column_list = list(df._get_numeric_data().columns)
    column_list = list(df.columns)
    categoric_columns = [i for i in column_list if
                         i not in numeric_column_list]
    # exclude random effect cols and SampleID
    # TODO do I want this in model? Fixed vs random effects?
    to_drop = [random_effect, sample_ID, "SampleID"]
    encode_this_list = [i for i in categoric_columns if i not in to_drop]

    x, z, clusters, y, feature_list, dummy_dict = \
        setup_df_do_encoding(df=df, random_effect=random_effect,
                             vars_to_encode=encode_this_list,
                             response_var=response_var, delta_flag=delta_flag,
                             join_flag=join_flag)

    if random_forest_type == "mixed":
        mrf, training_stats_plot = train_merf_model(x, z, clusters, y,
                                                    max_iters)
        y_predicted = mrf.predict(x, z, clusters)

        print(mean_squared_error(y, y_predicted))
        forest = mrf.trained_fe_model
        forest_score = round(forest.oob_score_, 3)
        print("Forest OOB score")
        print(forest_score)
        corr_plot = scatter_plot(y_predicted, y, forest_score)

        plot_list, feature_importance, top_features_list = shap_explainer(
            forest, x)
        print("\nranking of top 10 features\n")
        print(feature_importance.head(10))
        feature_importance.to_csv(out_file_prefix + "-feature-imp.txt",
                                  sep='\t', mode='a')
    else:
        # TODO
        mrf, training_stats_plot = train_rf_model(x, y)

    feature_selector = BorutaShap(model=forest,
                                  importance_measure='shap',
                                  classification=False)

    feature_selector.fit(X=x, y=y, n_trials=20, train_or_test="train",
                         sample=False)
    feature_selector.plot(which_features='all', figsize=(14, 10),
                          display=True)

    # find prefix for encoded values "ENC_Location_is_1_3" decoded is Location
    decoded_top_features_list = []
    for i in top_features_list:
        i_for_dict = i.split("_is_")[0] + "_is"
        try:
            decoded = dummy_dict[i_for_dict]
        except KeyError:
            decoded = i
        decoded_top_features_list.append(''.join(decoded))

    df = pd.DataFrame({'important.features': top_features_list,
                       'decoded.features': decoded_top_features_list})
    df.to_csv(out_file_prefix + "-top-features.txt", index=False, sep="\t")
    plot_partial_dependence(forest, x, features=top_features_list)
    plot_some_partial_dependence = plt.gcf()
    plot_some_partial_dependence.set_size_inches(width, height)
    plot_list.append(plot_some_partial_dependence)
    plot_list.append(training_stats_plot)
    plot_list.append(corr_plot)
    build_result_pdf(out_file, plot_list)


main(in_file=in_file, out_file=out_file, random_forest_type=random_forest_type,
     random_effect=random_effect, sample_ID=sample_ID,
     response_var=response_var, delta_flag=delta_flag, join_flag=join_flag)
