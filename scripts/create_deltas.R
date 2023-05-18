
#devtools::install_github("jbisanz/qiime2R")
# library(qiime2R)

build_datatable_viz <- NULL
source("scripts/install.R")
source("scripts/viz-datatable.R")

package_list <- c("dplyr", "usedist", "utils")
install_r_packages(package_list = package_list)

library(dplyr)
library(usedist)
library(utils)


# from config
output_folder <- snakemake@config[["out"]]

sample_id <- snakemake@config[["sample_id"]]

# change back to these names
input_timepoint <- snakemake@config[["timepoint"]]
input_study_id <- snakemake@config[["random_effect"]]


# from input, output, and params
in_file <- snakemake@input[["in_file"]]
out_file <- snakemake@output[["out_file"]]

reference_time <- snakemake@params[["reference_time"]]
absolute_values <- snakemake@params[["absolute_values"]]
build_datatable <- snakemake@params[["build_datatable"]]

# get name of distance matrix in params for rule, then
distance_matrices <- eval(parse(text = snakemake@params[["distance_matrices"]]))
distance_matrices <- lapply(distance_matrices, read.table, fill = TRUE)


include_distance_matrices <- function(distance_matrices, ref_sample, # nolint
                                      current_sample, df_new_comparison) {
  # add all distances (comparisons between sample timepoints); NA if none
    for (i in names(distance_matrices)){
      result <- tryCatch({
        result <-
          dist_subset(distance_matrices[[i]],
                      c(ref_sample[[sample_id]],
                        current_sample[[sample_id]]))[1]
      }, error = function(err) {
        result <- NA
        return(result)
      }, finally = function() {
        return(result)
      })
      df_new_comparison <- df_new_comparison %>%
        mutate({{i}} := result)
    }
  return(df_new_comparison) # update the new comparison
}



diffs_for_all_vars_per_subject <- function(delta_df, vars, sid_df,
                                           comparison, time, rt) {

  for (var in vars) {
    new_value <- sid_df[[var]][sid_df$timepoint_w_ == time]
    old_value <- sid_df[[var]][sid_df$timepoint_w_ == rt]
    if (typeof(sid_df[[var]]) != "character") {
      if (reference_time == "pairwise" && absolute_values == "yes") {
        # TODO convert to abs values at end? reduce computational time
        delta_df[[var]][delta_df$sid_delta == comparison] <-
        abs(new_value - old_value)
        delta_df[[paste0(var,
        "_reference")]][delta_df$sid_delta == comparison] <- old_value
      } else {
        delta_df[[var]][delta_df$sid_delta == comparison] <-
        new_value - old_value
        delta_df[[paste0(var,
        "_reference")]][delta_df$sid_delta == comparison] <- old_value
      }
    }
    if (typeof(sid_df[[var]]) == "character" || var == "timepoint_w_") {
      delta_df[[var]][delta_df$sid_delta == comparison] <-
        paste0(old_value, "__", new_value)
      delta_df[[paste0(var,
      "_reference")]][delta_df$sid_delta == comparison] <- paste0(old_value)
    } else {
      next
    }
  }

  return(delta_df)
}


all_reference_time_comparisons <- function(df, studyid, ref_times, time, vars, # nolint
                                           delta_df, dm_flag) {
  "timepoint_w_" <- "study_id_w_" <- "sid_delta" <- NULL
  for (rt in ref_times) {
    sid_df <- df %>%
      filter((timepoint_w_ == rt | timepoint_w_ == time) &
      study_id_w_ == studyid)
    # skip if no samples for subject at time points of interest
    # also skip if reference time is larger than time (would be duplicates)
    if (nrow(sid_df) < 2 || rt > time) {
      next
    }
    ref_sample <- sid_df %>%
      filter(study_id_w_ == {{studyid}} & timepoint_w_ == {{rt}})
    current_sample <- sid_df %>%
      filter(study_id_w_ == {{studyid}} & timepoint_w_ == {{time}})
    comparison <- paste0(studyid, "__", rt, "__", time)

    # Create new dataframe for each comparison
    df_new_comparison <- data.frame(sid_delta = comparison,
                                    timepoint_w_ = paste0(rt, "__", time),
                                    study_id_w_ = studyid)

    # if there are dms, add differences here for each subject
    if (dm_flag) {
      df_new_comparison <-
        include_distance_matrices(distance_matrices = distance_matrices,
                                  ref_sample = ref_sample,
                                  current_sample = current_sample,
                                  df_new_comparison = df_new_comparison)
    }

    delta_df <- dplyr::bind_rows(df_new_comparison, delta_df)

    # calculate differences for all variables for one subject
    delta_df <-
      diffs_for_all_vars_per_subject(delta_df = delta_df, vars = vars,
                                     sid_df = sid_df,
                                     comparison = comparison,
                                     time = time, rt = rt)
  }
  return(delta_df)
}


main <- function() {
  "timepoint_w_" <- "study_id_w_" <- "sid_delta" <- NULL
  dm_flag <- FALSE
  if (length(distance_matrices) > 0) {
    dm_flag <- TRUE
  }

  df <- read.csv(in_file, sep = "\t", check.names = FALSE)

  # FOR NSHAP study
  # df <- df %>% mutate(across(everything(), as.character))
  # df <- df %>% mutate(across(c(weight_sel, weight_adj, stratum,
  # drugs_count, np_count, cluster, gender, age, timepoint, happy), as.numeric))


  # if "timepoint_w_" or "study_id_w_" exist here, then can't go on
  if ("timepoint_w_" %in% colnames(df)) {
    print("Item is present in the List.")
  }
  if ("study_id_w_" %in% colnames(df)) {
    print("Item is present in the List.")
  }

  colnames(df)[colnames(df) == input_timepoint] <- "timepoint_w_"
  colnames(df)[colnames(df) == input_study_id] <- "study_id_w_"

  df$timepoint_w_ <- factor(df$timepoint_w_, ordered = TRUE)
  df$sid_delta <- NULL
  # get all variables/features, except ones used for script
  vars <- colnames(df %>% select(!study_id_w_))
  times <- unique(df$timepoint_w_)

  delta_df <- data.frame(sid_delta = character())
  # TODO ref_times are zero in 'previous' improve this.
  # TODO make function for ref_times

  # TODO add ability to use .QZA or TSV
  #unweighted.unifrac.file <- "data/unweighted_unifrac_distance_matrix.qza"
  #dm <- read_qza(unweighted.unifrac.file)
  #dm <- dm$data

  # calculate all reference time comparisons, per var, per person

  for (studyid in unique(df$study_id_w_)){
    for (time in times) {
      time <- as.numeric(time)
      # get different reference times to create deltas
      if (reference_time == "first") {
        ref_times <- 1
      } else if (reference_time == "previous") {
        ref_times <- time - 1
      } else if (reference_time == "pairwise") {
        # TODO times or timepoint
        ref_times <- as.numeric(levels(df$timepoint_w_))
        ref_times <- ref_times[ref_times != time]
      }
      delta_df <-
      all_reference_time_comparisons(df = df, studyid = studyid,
                                     ref_times = ref_times, time = time,
                                     vars = vars, delta_df = delta_df,
                                     dm_flag = dm_flag)
    }
  }

  delta_df <- delta_df %>%
    arrange(study_id_w_, timepoint_w_) %>%
    select(!sid_delta)

  colnames(delta_df)[colnames(delta_df) == "timepoint_w_"] <- input_timepoint
  colnames(delta_df)[colnames(delta_df) == "study_id_w_"] <- input_study_id

  delta_df <- delta_df %>%
    select(-ends_with("_reference"))

  write.table(delta_df, out_file, row.names = FALSE, sep = "\t")

  # TODO optional create visualizer
  # TODO change this to either df or load the dataframe
  if (build_datatable %in% c("TRUE", "True")) {
    input_file_name <-
    paste0(output_folder, "04-SELECTED-FEATURES-", reference_time,
           "/", reference_time, ".txt")
    output_file_name <- paste0(output_folder, "04-SELECTED-FEATURES-",
                               reference_time, "/vizualizer-", reference_time,
                               ".html")
    build_datatable_viz(input_file_name, output_file_name)
  }
}


main()