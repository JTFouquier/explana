


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
        in_file = "data/simulated-plants-fertilizer2.txt"
    output:
        out_file = "deltas/deltas-fertilizer-{reference}.txt"
    params:
        reference_time = "{reference}",
        absolute_values = "no",
    script:
        "scripts/create_deltas.R"

# rule integrate_data:
#     input:
#         in_file = .biom, .biom, .biom, etc
#         metadata = metadata file
#
#     output:
#         out_file = "merged-data.txt"

# Python in a complicated environment
# MERF has extremely specific Python requirements that could pose problems
# when combined with other software? Unsure.
# TODO add 'previous' or 'first' and other options to log file/report
# TODO prevent them from running pairwise from 'previous' file?
rule run_random_forest:
    input:
        in_file = "deltas/deltas-fertilizer-{reference}.txt"  # TODO above
    output:
        out_file = "random-forest/MERF-iter-20-deltas-{reference}.pdf"
    params:
        mixed_effects = "yes",
        random_effect = "StudyID",
        sample_ID = "StudyID.Timepoint",
        response_var = "ResponseMint",
        delta_flag = "yes"
    script:
        "scripts/random_forest.py"


rule run_random_forest_raw:
    input:
        in_file = "data/simulated-plants-fertilizer2.txt"  # TODO above
    output:
        out_file = "random-forest/MERF-iter-20-raw.pdf"
    params:
        mixed_effects = "yes",
        random_effect = "StudyID",
        sample_ID = "StudyID.Timepoint",
        response_var = "ResponseMint",
        delta_flag = "no"
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