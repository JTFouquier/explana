

import pandas as pd


def _dummy_var_creation(categorical_list, df, join_flag):
    # TODO SUBSET FLAG
    for i in categorical_list:
        subset_df = df[[i]]
        dum_df = pd.get_dummies(subset_df, columns=[i], prefix=["ENC_" + i +
                                                                "_is"])

        if join_flag:
            df = df.join(dum_df)

    if join_flag:
        return df
    else:
        return dum_df


def _dummy_this_list(categorical_list, df):
    for i in categorical_list:
        subset_df = df[[i]]
        dum_df = pd.get_dummies(subset_df, columns=[i], prefix=["ENC_" + i +
                                                                "_is"])

    return dum_df
