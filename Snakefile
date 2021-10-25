


# rule run_scnic:
# Python in its own environment
#
#     input:
#         script =
#         data =

# look at results
# rule scnic_between:
# # ideas to test
# R simple environment
# this step does not have to be done
# TODO add appended name to everything idea
# TODO action items
rule create_deltas:
    input:
        in_file = "data/simulated-diet.txt"
    output:
        out_file = "deltas/simulated-diet-deltas-{reference}.txt"
    params:
        reference_time = "{reference}",
        absolute_values = "no",
    script:
        "scripts/create_deltas.R"

# rule integrate_data:
#     input:
#         in_file = .biom, .txt, (and multiples of these) etc
#         metadata = metadata file
#
#     output:
#         out_file = "merged-data.txt"

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
rule run_random_forest:
    input:
        in_file = "deltas/simulated-diet-deltas-{reference}.txt"  # TODO above
    output:
        out_file = "random-forest/{mixed}-RF-{deltas}-{re_timepoint}-{reference}.pdf"
    params:
        random_forest_type = "{mixed}", # mixed or fixed
        random_effect = "StudyID",
        sample_ID = "StudyID.Timepoint",
        response_var = "IL6",
        delta_flag = "{deltas}", # raw or deltas
        iterations = 8, # iterations, 20 is suggested, 10 for testing
        re_timepoint = "{re_timepoint}" # re_timepoint or no_re
    script:
        "scripts/random_forest.py"

# run random forest on repeated studies data, no deltas, still mixed effects
rule run_random_forest_raw:
    input:
        in_file = "data/simulated-diet.txt"  # TODO above
    output:
        out_file = "random-forest/MERF-iter-8-mixed-raw.pdf"
    params:
        random_forest_type = "mixed",
        random_effect = "StudyID",
        sample_ID = "StudyID.Timepoint",
        response_var = "IL6",
        delta_flag = "raw",
        iterations = 8, # iterations, 20 is suggested, 10 for testing
        re_timepoint = "no_re" # re_timepoint or no_re
    script:
        "scripts/random_forest.py"


rule run_post_hoc_stats:
    input:
        deltas_first = "deltas/simulated-diet-deltas-first.txt",
        deltas_previous = "deltas/simulated-diet-deltas-previous.txt",
        deltas_pairwise = "deltas/simulated-diet-deltas-pairwise.txt",
        df_raw = "data/simulated-diet.txt",
        fixed_effects_first = "random-forest/mixed-RF-deltas-re_timepoint-first-top-features.txt",
        fixed_effects_previous = "random-forest/mixed-RF-deltas-re_timepoint-previous-top-features.txt",
        fixed_effects_pairwise = "random-forest/mixed-RF-deltas-re_timepoint-pairwise-top-features.txt",
        fixed_effects_raw = "random-forest/MERF-iter-8-mixed-raw-top-features.txt"
    output:
        out_file = "post-hoc-analysis.html"
    params:
        response_var = "IL6"
    script:
        "scripts/post_hoc_tests.R"

# look at results
#
#
#
# find correlating vars after TODO
#
# MERF
# Python 3.6 - 3.7