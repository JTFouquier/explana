

configfile: "config.yaml"

rule all:
    input:
        first = config["out"] + "04-SELECTED-FEATURES-first/mixed-RF-deltas-"
                "re_timepoint-first.pdf",
        previous = config["out"] + "04-SELECTED-FEATURES-previous/mixed-RF-deltas-"
                   "re_timepoint-previous.pdf",
        pairwise = config["out"] + "04-SELECTED-FEATURES-pairwise/mixed-RF-deltas-"
                   "re_timepoint-pairwise.pdf",
        original = config["out"] + "04-SELECTED-FEATURES-original/MERF-original.pdf",


rule dim_reduction_pca:
    input:
        in_file = "data/real-data-no-asvs.txt",
    output:
        out_file = config["out"] + "01-DIM-PCA-metadata/"
                   "final-pca-dim-reduction.txt"
    params:
        pca_groups_list = "list(list("
                          "pca_name='PCA1',"
                          "feature_list=c('Triglycerides', 'LDL', 'Leptin', "
                          "'Adiponectin', 'Insulin', 'Glucose'),"
                          "cum_var_thres=70), list(pca_name='PCA2',"
                          "feature_list=c('Bact', 'Prev'), "
                          "cum_var_thres=70))"
    script:
        "scripts/dim_reduction_pca.R"

# TODO remove in_file_metadata
rule dim_reduction_scnic:
    input:
        in_file = "data/real-feature-table.biom",
        in_file_metadata = "data/real-data-no-asvs.txt"
    output:
        out_file = config["out"] + "01-DIM-SCNIC-biom/"
                   "SCNIC_modules.txt"
    script:
        "scripts/dim_reduction_scnic.py"

# TODO have the default be to merge all datasets output from 01 steps
# TODO add drop columns or rows here instead
rule integrate_datasets:
    input:
        dataset_list = [config["out"] + "01-DIM-PCA-metadata/"
                        "final-pca-dim-reduction.txt",
                        config["out"] + "01-DIM-SCNIC-biom/"
                        "SCNIC_modules.txt"]
    output:
        out_file = config["out"] + "02-MERGED-DATA/"
                   "final-merged-dfs.txt"
    script:
        "scripts/merge_datasets.py"


rule make_delta_datasets:
    input:
        in_file = config["out"] + "02-MERGED-DATA/final-merged-dfs.txt"
    output:
        out_file = config["out"] + "04-SELECTED-FEATURES-{reference}/deltas-{reference}.txt"
    params:
        reference_time = "{reference}",
        absolute_values = "no",
        build_visualizer = True, # TODO have this be optional default true
    script:
        "scripts/create_deltas.R"


rule visualize_datasets:
    input:
        in_file = config["out"] + "04-SELECTED-FEATURES-{reference}/"
                                  "deltas-{reference}.txt"
    output:
        out_file = config["out"] + "04-SELECTED-FEATURES-{reference}/"
                                   "vizualizer-{reference}.html"
    script:
        "scripts/viz-datatable.R"


# Python in a complicated environment
# MERF has extremely specific Python requirements that could pose problems
# when combined with other software? Unsure.
# TODO add 'previous' or 'first' and other options to log file/report
# TODO prevent them from running pairwise from 'previous' file?
# TODO fix underscores in names.
#
# Also, this whole thing is done a little strange with the variables in the
# name of the output. Not sure how I want to handle this, but it was because
# I wanted my file names to reflect the arguments

# TODO check constrain, and wildcards (not a good use)
rule random_forest_deltas:
    input:
        in_file = config["out"] + "04-SELECTED-FEATURES-{reference}/deltas-{reference}.txt"
    output:
        out_file = config["out"] + "04-SELECTED-FEATURES-{reference}/{mixed}-RF-{deltas}-"
                   "{re_timepoint}-{reference}.pdf"
    params:
        random_forest_type = "{mixed}",  # mixed or fixed
        random_effect = config["random_effect"],
        sample_ID = "StudyID.Timepoint",
        drop_rows = {"SexualClassification": "Women_Women"},
        constrain_rows = {"Diet": "Agrarian_Agrarian"},
        drop_cols = ["Inflammation","TotalCholesterol", "PCA", "HOMAIR"],
        constrain_cols = [],
        response_var = config["response_var"],
        delta_flag = "{deltas}",  # raw or deltas
        iterations = 5,  # iterations, 20 is suggested, 10 for testing
        re_timepoint = "{re_timepoint}"  # re_timepoint or no_re
    script:
        "scripts/random_forest.py"


rule random_forest_original:
    input:
        in_file = config["out"] + "02-MERGED-DATA/final-merged-dfs.txt"
    output:
        # report("{my_dir}/04-SELECTED-FEATURES-original/MERF-original.pdf",
        #         caption="report/plot_pcas.rst", category="Selected Features",
        #         subcategory="Original Dataset"),
        out_file = config["out"] + "04-SELECTED-FEATURES-original/MERF-original.pdf"
    params:
        random_forest_type = "mixed",
        random_effect = config["response_var"],
        sample_ID = "StudyID.Timepoint",
        drop_rows = {"SexualClassification": "Women"},
        constrain_rows = {"Diet": "Agrarian"},
        drop_cols = ["Inflammation", "TotalCholesterol", "PCA", "HOMAIR"],
        constrain_cols = [],
        response_var = config["response_var"],
        delta_flag = "raw",
        iterations = 5,
        re_timepoint = "no_re"
    script:
        "scripts/random_forest.py"

# rule TODO might be better to build report dynamically:
#     output:
#         report(
#             directory("random-forest/TEST/HDL/01-DIM-PCA-metadata/PCA1/"),
#             patterns=["{name}.pdf"],
#             caption="report/plot_pcas.rst")


# TODO program if statements here and make report reflect missing data

rule render_report:
    input:
        original = rules.random_forest_original.output.out_file,
        original_shap= config["out"] + "04-SELECTED-FEATURES-original/MERF-original-"
                       "accepted_SHAP.svg",
        original_boruta = config["out"] + "04-SELECTED-FEATURES-original/MERF-original-boruta-"
                          "accepted-features.svg",
        original_log = config["out"] + "04-SELECTED-FEATURES-original/MERF-original-log.txt",

        first = config["out"] + "04-SELECTED-FEATURES-first/mixed-RF-deltas-"
                "re_timepoint-first.pdf",
        first_shap= config["out"] + "04-SELECTED-FEATURES-first/mixed-RF-deltas-re_timepoint-"
                    "first-accepted_SHAP.svg",
        first_boruta= config["out"] + "04-SELECTED-FEATURES-first/mixed-RF-deltas-re_timepoint-"
                      "first-boruta-accepted-features.svg",
        first_log = config["out"] + "04-SELECTED-FEATURES-first/mixed-RF-deltas-"
                    "re_timepoint-first-log.txt",

        previous = config["out"] + "04-SELECTED-FEATURES-previous/mixed-RF-deltas-"
                   "re_timepoint-previous.pdf",
        previous_shap= config["out"] + "04-SELECTED-FEATURES-previous/mixed-RF-deltas-"
                       "re_timepoint-previous-accepted_SHAP.svg",
        previous_boruta = config["out"] + "04-SELECTED-FEATURES-previous/mixed-RF-deltas-"
                          "re_timepoint-previous-boruta-accepted-features.svg",
        previous_log= config["out"] + "04-SELECTED-FEATURES-previous/mixed-RF-deltas-"
                      "re_timepoint-previous-log.txt",

        pairwise = config["out"] + "04-SELECTED-FEATURES-pairwise/mixed-RF-deltas-"
                   "re_timepoint-pairwise.pdf",
        pairwise_shap= config["out"] + "04-SELECTED-FEATURES-pairwise/mixed-RF-deltas-"
                       "re_timepoint-pairwise-accepted_SHAP.svg",
        pairwise_boruta = config["out"] + "04-SELECTED-FEATURES-pairwise/mixed-RF-deltas-"
                          "re_timepoint-pairwise-boruta-accepted-features.svg",
        pairwise_log= config["out"] + "04-SELECTED-FEATURES-previous/mixed-RF-deltas-"
                      "re_timepoint-previous-log.txt",

        dag_plot = "scripts/dag.svg",
    output:
        md_doc="report.html",
    params:
        iters=rules.random_forest_original.params.iterations,
        response_var= config["response_var"],
        drop_cols= rules.random_forest_original.params.drop_cols,
        constrain_rows= rules.random_forest_original.params.constrain_rows,
        random_effect= rules.random_forest_original.params.random_effect,
        drop_rows= rules.random_forest_original.params.drop_rows,
    script:
        "scripts/report.Rmd"

    # TODO needs work
# rule run_post_hoc_stats:
#     pass
#     input:
#         deltas_first = "random-forest/TEST/HDL/04-feature-selection-HDL-"
#                        "agrarian-no-women-first/deltas-first.txt",
#         deltas_previous = "random-forest/TEST/HDL/04-feature-selection-HDL-"
#                           "agrarian-no-women-previous/deltas-previous.txt",
#         deltas_pairwise = "random-forest/TEST/HDL/04-feature-selection-HDL-"
#                           "agrarian-no-women-pairwise/deltas-pairwise.txt",
#         df_raw = "random-forest/TEST/HDL/02-data-integration/"
#                  "final-merged-dfs.txt",
#         fixed_effects_first = "random-forest/ptests/real-HDLAgrarian-no-women-no-asvs-first/mixed-RF-deltas-re_timepoint-first-top-features.txt",
#         fixed_effects_previous = "random-forest/ptests/real-HDLAgrarian-no-women-no-asvs-previous/mixed-RF-deltas-re_timepoint-previous-top-features.txt",
#         fixed_effects_pairwise = "random-forest/ptests/real-HDLAgrarian-no-women-no-asvs-pairwise/mixed-RF-deltas-re_timepoint-pairwise-top-features.txt",
#         fixed_effects_raw = "random-forest/ptests/real-HDLAgrarian-no-women-no-asvs-raw/MERF-iter-20-mixed-raw-top-features.txt"
#     output:
#         out_file = "random-forest/post-hoc-analysis-real-HDLAgrarian-no-women-no-asvs.html"
#     params:
#         response_var = "HDL"
#     script:
#         "scripts/post_hoc_tests.R"
