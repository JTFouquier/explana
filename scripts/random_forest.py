

import sys
import shap
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from merf.merf import MERF
# from merf.viz import plot_merf_training_stats
from sklearn.inspection import plot_partial_dependence

from hot_encoding import _dummy_var_creation, _dummy_this_list
from matplotlib.backends.backend_pdf import PdfPages

sns.set_context("talk")
from scipy import stats

from sklearn.ensemble import RandomForestRegressor

# reading
# https://slundberg.github.io/shap/notebooks/plots/dependence_plot.html
# https://shap.readthedocs.io/en/latest/example_notebooks/tabular_examples/tree_based_models/Force%20Plot%20Colors.html?highlight=force%20plot
# https://www.sciencedirect.com/science/article/pii/S1532046418301758

"""
TODO list 
Add baseline information and perform tests
"""

# access values from snakemake's Snakefile
in_file = snakemake.input["in_file"]
out_file = snakemake.output["out_file"]

mixed_effects = snakemake.params["mixed_effects"]
random_effect = snakemake.params["random_effect"]
sample_ID = snakemake.params["sample_ID"]
response_var = snakemake.params["response_var"]
delta_flag = snakemake.params["delta_flag"]

# set graphics dimensions
width = 15
height = 10
max_display = 8
top_features = 4

def setup_df_do_encoding(df_in, random_effect, encode, rv, delta_flag):
    df_in = _dummy_var_creation(encode, df_in, join_flag=True)
    # TODO fix hardcoding & add sample ID for all
    if delta_flag == "yes":
        drop_cols = ["StudyID.Timepoint"] + [random_effect] + encode + [rv]
    else:
        drop_cols = [random_effect] + encode + [rv]

    dum = df_in.drop(drop_cols, axis=1)
    feature_list = dum.columns

    df_in.to_csv("dummy-variables-df.txt", index=False, sep="\t")
    dum.to_csv("dummy-variables-drop.txt", index=False, sep="\t")

    x = df_in[dum.columns]
    # z = sunlight_numpy # TODO look at z matrix
    z = np.ones((len(x), 1))
    clusters_train = df_in["StudyID"]
    y = df_in[rv]

    return x, z, clusters_train, y, feature_list


def merf_training_stats_plot(mrf, num_clusters):
    # print training stats TODO understand this plot more
    # plot_merf_training_stats(mrf, num_clusters_to_plot=num_clusters)
    training_stats_plot = plt.gcf()
    plt.close()
    return training_stats_plot


def y_estimate_y_known_correlation(y_estimated, y):
    corr = stats.pearsonr(y_estimated, y)
    log_string = "Predicted vs observed response R2: " + \
                 str(round((corr[0] * corr[0]), 3)) + "p-value: " + \
                 str(round(corr[1], 3))
    return log_string


def train_MERF_model(x, z, clusters, y):
    # TODO CHANGE THIS later! NEEDS MORE ITERATIONS
    mrf = MERF(max_iterations=10)
    mrf.fit(x, z, clusters, y)
    training_stats_plot = merf_training_stats_plot(mrf, 7)

    return mrf, training_stats_plot


def train_RF_model(x, y):

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


def important_shapley_features(shap_values, X):
    feature_names = list(X.columns.values)
    vals = np.abs(shap_values).mean(0)
    feature_importance = pd.DataFrame(list(zip(feature_names, vals)),
                                      columns=['col_name',
                                               'feature_importance_vals'])
    feature_importance.sort_values(by=['feature_importance_vals'],
                                   ascending=False, inplace=True)
    print(feature_importance)


def shap_explainer(trained_forest_model, X):
    explainer = shap.TreeExplainer(trained_forest_model)
    shap_values = explainer.shap_values(X)

    plot_list = []

    shap.summary_plot(shap_values, X, show=False, plot_size=(width, height),
                      max_display=max_display)
    plot_list.append(graphic_handling(plt.gcf()))

    # Bar plot, just a different way to look
    shap.summary_plot(shap_values, X, plot_type="bar", show=False,
                      plot_size=(width, height), max_display=max_display)
    plot_list.append(graphic_handling(plt.gcf()))

    feature_names = list(X.columns.values)
    vals = np.abs(shap_values).mean(0)
    feature_importance = pd.DataFrame(list(zip(feature_names, vals)),
                                      columns=['col_name',
                                               'feature_importance_vals'])
    feature_importance.sort_values(by=['feature_importance_vals'],
                                   ascending=False, inplace=True)
    print(feature_importance)

    top_features_list = list(feature_importance['col_name'][0:top_features])
    for col in top_features_list:
        shap.dependence_plot(col, shap_values, X, show=False)
        # plot_list.append(graphic_handling(plt.gcf()))
        p = plt.gcf()
        p.set_size_inches(width, height)
        plot_list.append(p)
        plt.close()

    for col in top_features_list:
        i = X.columns.get_loc(col)
        shap.force_plot(explainer.expected_value, shap_values[i, :],
                        X.iloc[i, :], show=False, matplotlib=True)
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


def main(in_file, out_file, mixed_effects, random_effect, sample_ID,
         response_var, delta_flag):
    df = pd.read_csv(in_file, sep="\t")
    numeric_column_list = list(df._get_numeric_data().columns)
    column_list = list(df.columns)
    categoric_columns = [i for i in column_list if
                         i not in numeric_column_list]

    # exclude random effect cols and SampleID
    # TODO do I want this in model? Fixed vs random effects?
    to_drop = [random_effect, sample_ID]
    encode = [i for i in categoric_columns if i not in to_drop]

    X, Z, clusters, Y, feature_list = \
        setup_df_do_encoding(df_in=df, random_effect=random_effect,
                             encode=encode, rv=response_var,
                             delta_flag=delta_flag)

    if mixed_effects == "yes":
        mrf, training_stats_plot = train_MERF_model(X, Z, clusters, Y)
        y_predicted = mrf.predict(X, Z, clusters)
        mrf = mrf.trained_rf  # TODO why is this trained fixed effects
        # mrf = mrf.trained_fe_model # TODO this is old or new in MERF tool?

    else:
        mrf, training_stats_plot = train_RF_model(X, Y)
        y_predicted = mrf.predict(X)

    log_info = y_estimate_y_known_correlation(y_predicted, Y)
    f_out = open(out_file + "Log.txt", "w")
    f_out.write(log_info)
    f_out.close()

    plot_list, feature_importance, top_features_list = shap_explainer(mrf, X)

    feature_importance.to_csv(out_file + "feature_importance.txt",
                              sep='\t', mode='a')

    plot_partial_dependence(mrf, X, features=top_features_list)
    plot_some_partial_dependence = plt.gcf()
    plot_some_partial_dependence.set_size_inches(width, height)
    plt.close()
    plot_list.append(plot_some_partial_dependence)
    build_result_pdf(out_file, plot_list)


main(in_file=in_file, out_file=out_file, mixed_effects=mixed_effects,
     random_effect=random_effect, sample_ID=sample_ID,
     response_var=response_var, delta_flag=delta_flag)
