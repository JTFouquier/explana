
from matplotlib.backends.backend_pdf import PdfPages


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

def _build_result_pdf(out_file, out_file_prefix, plot_list,
                      plot_file_name_list):
    pp = PdfPages(out_file)
    my_count = 0
    # TODO add only important plot list (don't save svg or pdf)
    for i in plot_list:
        file_name = out_file_prefix + plot_file_name_list[my_count]
        my_count += 1

        if "dependence-plot" in file_name:
            pass
        else:
            i.savefig(file_name + ".svg", format='svg')

        pp.savefig(i)
    pp.close()
