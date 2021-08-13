


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
        in_file = "data/simulated-plants.txt"
    output:
        out_file = "deltas/simulated-plants-deltas-{reference}.txt"
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
        in_file = "deltas/simulated-plants-deltas-{reference}.txt"  # TODO above
    output:
        out_file = "random-forest/{mixed}-RF-{deltas}-{re_timepoint}-{reference}.pdf"
    params:
        random_forest_type = "{mixed}", # mixed or fixed
        random_effect = "StudyID",
        sample_ID = "StudyID.Timepoint",
        response_var = "ResponseMint",
        delta_flag = "{deltas}", # raw or deltas
        iterations = 20, # iterations, 20 is suggested, 10 for testing
        re_timepoint = "{re_timepoint}" # re_timepoint or no_re
    script:
        "scripts/random_forest.py"

# run random forest on repeated studies data, no deltas, still mixed effects
rule run_random_forest_raw:
    input:
        in_file = "data/simulated-plants.txt"  # TODO above
    output:
        out_file = "random-forest/MERF-iter-20-mixed-raw.pdf"
    params:
        random_forest_type = "mixed",
        random_effect = "StudyID",
        sample_ID = "StudyID.Timepoint",
        response_var = "ResponseMint",
        delta_flag = "raw",
        iterations = 20, # iterations, 20 is suggested, 10 for testing
        re_timepoint = "no_re" # re_timepoint or no_re
    script:
        "scripts/random_forest.py"

# look at results
#
#
#
# find correlating vars after TODO
#
# MERF
# Python 3.6 - 3.7