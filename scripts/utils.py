
from matplotlib.backends.backend_pdf import PdfPages


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
