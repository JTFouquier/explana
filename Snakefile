


# TODO add appended name to everything idea
# TODO action items
rule all:
    input:
        first = "random-forest/TEST/HDL/04-SELECTED-FEATURES--HDL-"
                "agrarian-no-women-first/mixed-RF-deltas-"
                "re_timepoint-first.pdf",
        previous = "random-forest/TEST/HDL/04-SELECTED-FEATURES--"
                   "HDL-agrarian-no-women-previous/mixed-RF-deltas-"
                   "re_timepoint-previous.pdf",
        pairwise = "random-forest/TEST/HDL/04-SELECTED-FEATURES--"
                   "HDL-agrarian-no-women-pairwise/mixed-RF-deltas-"
                   "re_timepoint-pairwise.pdf",
        original = "random-forest/TEST/HDL/04-SELECTED-FEATURES--"
                   "HDL-agrarian-no-women-original/MERF-original.pdf"
                   "HDL-agrarian-no-women-original/MERF-original.pdf",
        viz_first = "random-forest/TEST/HDL/04-SELECTED-FEATURES--HDL-"
                    "agrarian-no-women-first/vizualizer-first.html",
        viz_previous = "random-forest/TEST/HDL/04-SELECTED-FEATURES--HDL-"
                       "agrarian-no-women-previous/vizualizer-previous.html",
        viz_pairwise = "random-forest/TEST/HDL/04-SELECTED-FEATURES--HDL-"
                       "agrarian-no-women-pairwise/vizualizer-pairwise.html"


rule dim_reduction_pca:
    input:
        in_file = "data/real-data-no-asvs.txt",
    output:
        out_file = "random-forest/TEST/HDL/01-DIM-PCA--metadata/"
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
        out_file = "random-forest/TEST/HDL/01-DIM-SCNIC--biom/"
                   "SCNIC_modules.txt"
    script:
        "scripts/dim_reduction_scnic.py"

# TODO have the default be to merge all datasets output from 01 steps
rule merge_datasets:
    input:
        dataset_list = ["random-forest/TEST/HDL/01-DIM-PCA--metadata/"
                        "final-pca-dim-reduction.txt",
                        "random-forest/TEST/HDL/01-DIM-SCNIC--biom/"
                        "SCNIC_modules.txt"]
    output:
        out_file = "random-forest/TEST/HDL/02-MERGED-DATA/"
                   "final-merged-dfs.txt"
    script:
        "scripts/merge_datasets.py"


rule make_delta_datasets:
    input:
        in_file = "random-forest/TEST/HDL/02-MERGED-DATA/"
                  "final-merged-dfs.txt"
    output:
        out_file = "random-forest/TEST/HDL/04-SELECTED-FEATURES--HDL-"
                   "agrarian-no-women-{reference}/deltas-{reference}.txt"
    params:
        reference_time = "{reference}",
        absolute_values = "no",
    script:
        "scripts/create_deltas.R"


rule visualize_datasets:
    input:
        in_file = "random-forest/TEST/HDL/04-SELECTED-FEATURES--HDL-"
                  "agrarian-no-women-{reference}/deltas-{reference}.txt"
    output:
        out_file = "random-forest/TEST/HDL/04-SELECTED-FEATURES--HDL-"
                   "agrarian-no-women-{reference}/vizualizer-{reference}.html"
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
        in_file = "random-forest/TEST/HDL/04-SELECTED-FEATURES--HDL-"
                  "agrarian-no-women-{reference}/deltas-{reference}.txt"
    output:
        out_file = "random-forest/TEST/HDL/04-SELECTED-FEATURES--HDL-"
                   "agrarian-no-women-{reference}/{mixed}-RF-{deltas}-"
                   "{re_timepoint}-{reference}.pdf"
    params:
        random_forest_type = "{mixed}",  # mixed or fixed
        random_effect = "StudyID",
        sample_ID = "StudyID.Timepoint",
        drop_rows = {"SexualClassification": "Women_Women"},
        constrain_rows = {"Diet": "Agrarian_Agrarian"},
        drop_cols = ["Inflammation","TotalCholesterol", "PCA", "HOMAIR"],
        constrain_cols = [],
        response_var = "HDL",
        delta_flag = "{deltas}",  # raw or deltas
        iterations = 20,  # iterations, 20 is suggested, 10 for testing
        re_timepoint = "{re_timepoint}"  # re_timepoint or no_re
    script:
        "scripts/random_forest.py"


rule random_forest_original:
    input:
        in_file = "random-forest/TEST/HDL/02-MERGED-DATA/"
                  "final-merged-dfs.txt"
    output:
        out_file = "random-forest/TEST/HDL/04-SELECTED-FEATURES--HDL-"
                   "agrarian-no-women-original/MERF-original.pdf"
    params:
        random_forest_type = "mixed",
        random_effect = "StudyID",
        sample_ID = "StudyID.Timepoint",
        drop_rows = {"SexualClassification": "Women"},
        constrain_rows = {"Diet": "Agrarian"},
        drop_cols = ["Inflammation", "TotalCholesterol", "PCA", "HOMAIR"],
        constrain_cols = [],
        response_var = "HDL",
        delta_flag = "raw",
        iterations = 20,
        re_timepoint = "no_re"
    script:
        "scripts/random_forest.py"

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

# TODO make report
# rule make_report:
#     input:
#         "path/to/inputfile",
#         "path/to/other/inputfile"
#     output:
#         "path/to/report.html",
#     script:
#         "path/to/report.Rmd"