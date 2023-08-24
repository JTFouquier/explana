


"""
Adapted from:
https://github.com/Ekeany/Boruta-Shap/blob/master/src/BorutaShap.py

Modified by Jennifer Fouquier to remove plt.show() or plt.close() so I can
save to PDF file report and customize figure titles for improved report.


"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import re


def _create_mapping_of_features_to_attribute(self, maps=[]):

    rejected = list(self.rejected)
    tentative = list(self.tentative)
    accepted = list(self.accepted)
    shadow = ['Max_Shadow', 'Median_Shadow', 'Min_Shadow', 'Mean_Shadow']

    tentative_map = self.create_list(tentative, maps[0])
    rejected_map = self.create_list(rejected, maps[1])
    accepted_map = self.create_list(accepted, maps[2])
    shadow_map = self.create_list(shadow, maps[3])

    values = tentative_map + rejected_map + accepted_map + shadow_map
    keys = tentative + rejected + accepted + shadow

    return self.to_dictionary(keys, values)


def _make_box_plot(self, data, X_rotation, X_size, y_scale, figsize,
                   which_features):

    if y_scale == 'log':
        minimum = data['value'].min()
        if minimum <= 0:
            data['value'] += abs(minimum) + 0.01

    order = data.groupby(by=["Methods"])["value"].mean().sort_values(
        ascending=False).index
    my_palette = \
        _create_mapping_of_features_to_attribute(self,
                                                 maps=['yellow', 'red',
                                                       'green', 'blue'])
    # Use a color palette
    # dynamically change width of figure if there are many values
    total_x_vals = len(set(data["Methods"]))
    allowed_x_per_width = 80
    x_ratio = total_x_vals/allowed_x_per_width
    if x_ratio < 1:
        boruta_width = figsize[0]
    else:
        boruta_width = figsize[0]*x_ratio
    figsize = (boruta_width, figsize[1])
    plt.figure(figsize=figsize)
    ax = sns.boxplot(x=data["Methods"], y=data["value"], order=order,
                     palette=my_palette)

    if y_scale == 'log':
        ax.set(yscale="log")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=X_rotation, size=X_size)

    for index, label in enumerate(ax.xaxis.get_ticklabels()):
        new_label = re.split(r", |\)", str(label))[2]
        new_label = new_label.replace("'", "")

    ax.set_title(which_features.capitalize() + " Features from BorutaShap",
                 fontsize=12)
    ax.set_ylabel('Z-Score', fontsize=10)
    ax.set_xlabel('Features', fontsize=10)
    return boruta_width


def _boruta_shap_plot(self, X_rotation=90, X_size=9, figsize=(13, 8),
                      y_scale='log', which_features='all'):
    """
    creates a boxplot of the feature importances
    Parameters
    ----------
    X_rotation: int
        Controls the orientation angle of the tick labels on the X-axis
    X_size: int
        Controls the font size of the tick labels
    y_scale: string
        Log transform of the y axis scale as hard to see the plot as it is
        normally dominated by two or three
        features.
    which_features: string
        Despite efforts if the number of columns is large the plot becomes
        cluttered so this parameter allows you to
        select subsets of the features like the accepted, rejected or tentative
        features default is all.
        controls if the output is displayed or not, set to false when running
        test scripts
    """
    # data from wide to long
    data = self.history_x.iloc[1:]
    data['index'] = data.index
    data = pd.melt(data, id_vars='index', var_name='Methods')

    decision_mapper = self.create_mapping_of_features_to_attribute \
        (maps=['Tentative', 'Rejected', 'Accepted', 'Shadow'])
    data['Decision'] = data['Methods'].map(decision_mapper)
    data.drop(['index'], axis=1, inplace=True)

    options = {'accepted': self.filter_data(data, 'Decision', 'Accepted'),
               'tentative': self.filter_data(data, 'Decision', 'Tentative'),
               'rejected': self.filter_data(data, 'Decision', 'Rejected'),
               'all': data}

    self.check_if_which_features_is_correct(which_features)
    data = options[which_features.lower()]

    boruta_width = _make_box_plot(self, data=data, X_rotation=X_rotation,
                                  X_size=X_size, y_scale=y_scale,
                                  figsize=figsize,
                                  which_features=which_features)
    return boruta_width

