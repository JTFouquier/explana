import pandas as pd
import re
import urllib.parse

feature_file = snakemake.input["feature_file"]
response = snakemake.config["response_var"]
out_file = snakemake.output["out_file"]


def perform_pubmed_search(search_terms):
    encoded_search_terms = urllib.parse.quote(search_terms)
    base_url = 'https://pubmed.ncbi.nlm.nih.gov/?term='
    url = base_url + encoded_search_terms
    return (url)


def build_url_list(feature_file, out_file, compiled_flag):
    # summary needs compiled url file with all datasets
    # not compiled is for individual datasets
    if compiled_flag:
        col_names = ["dataset", "important_features"]
    if not compiled_flag:
        col_names = ["important_features"]
    try:
        df = pd.read_csv(feature_file, sep="\t", usecols=col_names)
    except:
        with open(out_file, "w") as file:
            file.write("Failed Analysis")
        return

    def make_url(x):
        # if people use numbers in column names convert to string; fix earlier
        x = str(x)
        try:
            feature = re.split("_is_", x)[0]
            feature = re.split("ENC_", feature)[1]

            feature_values = re.split("_is_", x)[1]
            try:
                feature_value_1 = re.split("__", feature_values)[0]
                feature_value_2 = re.split("__", feature_values)[1]
                url = perform_pubmed_search(feature_value_1 + " AND " +
                                            feature_value_2 + " AND " +
                                            feature + " AND " + response)
            except:
                url = perform_pubmed_search(feature_values + " AND " + feature
                                            + " AND " + response)
        except IndexError:
            url = perform_pubmed_search(x + " AND " + response)

        return (url)

    df_features = df["important_features"]
    df_url = df_features.apply(make_url).to_frame()
    df_url = df_url.rename(columns={"important_features": 'url'})

    # update df with other cols
    df_final = pd.concat([df, df_url.reindex(df.index)], axis=1)

    df_final.to_csv(out_file, sep="\t", index=False)


def main(feature_file, out_file):

    out_p = str(snakemake.config["out"])

    build_url_list(feature_file=feature_file, out_file=out_file,
                   compiled_flag=True)

    ref_list = ["original", "first", "previous", "pairwise"]
    for i in range(0, len(ref_list)):
        i = ref_list[i]
        ref_features = snakemake.input[i + "_features"]

        out_file = out_p + str(snakemake.config["path_" + i]) + i + "-urls.txt"
        build_url_list(feature_file=ref_features, out_file=out_file,
                       compiled_flag=False)


main(feature_file=feature_file, out_file=out_file)
