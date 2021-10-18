

import pandas as pd


def _dummy_var_creation(categorical_list, df, join_flag):
    # TODO SUBSET FLAG
    dummy_dict = {}
    for orig_name in categorical_list:
        subset_df = df[[orig_name]]
        encoded_name = "ENC_" + orig_name + "_is"
        dum_df = pd.get_dummies(subset_df, columns=[orig_name],
                                prefix=[encoded_name])
        dummy_dict[encoded_name] = [orig_name]
        # Keep adding encoded vars to df before returning
        if join_flag:
            df = df.join(dum_df)

    if join_flag:
        return df, dummy_dict

    # TODO without join flag... revisit this idea
    else:
        return dum_df