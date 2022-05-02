

def _subset_simple(df, drop_rows, subset_rows, drop_cols, subset_cols):

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

    return df