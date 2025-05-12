
source("scripts/install.R")
package_list <- c("usedist", "utils", "tidyverse")
install_r_packages(package_list = package_list)

library(tidyverse)
library(usedist)
library(utils)


# from config
output_folder <- snakemake@config[["out"]]
response_var <- snakemake@config[["response_var"]]

sample_id <- snakemake@config[["sample_id"]]

# change back to these names
input_timepoint <- "timepoint_explana"
input_study_id <- snakemake@config[["random_effect"]]

analyze_first <- snakemake@config[["analyze_first"]]
analyze_previous <- snakemake@config[["analyze_previous"]]
analyze_pairwise <- snakemake@config[["analyze_pairwise"]]
include_reference_values <- snakemake@params[["include_reference_values"]]
absolute_values <- snakemake@params[["absolute_values"]]

# from input, output, and params
in_file <- snakemake@input[["in_file"]]
out_file <- snakemake@output[["out_file"]]

reference_time <- snakemake@params[["reference_time"]]

# get name of distance matrix in params for rule, then load by name
distance_matrices <- eval(parse(text = snakemake@params[["distance_matrices"]]))

open_dms <- function(dm) {

  dm <- read.table(dm, sep = "\t", header = TRUE, check.names = FALSE,
  row.names = 1)
  # dm is already in distance format, so load as matrix
  dm <- as.matrix(dm, dimnames = colnames(df))
  return(dm)
}

distance_matrices <- lapply(distance_matrices, open_dms)


include_distance_matrices <- function(distance_matrices, ref_sample, # nolint
                                      current_sample, df_new_comparison) {
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
        dplyr::mutate({{i}} := result)
    }
  return(df_new_comparison) # update the new comparison
}


diffs_for_all_vars_per_subject <- function(delta_df, vars, sid_df,
                                           comparison, time, rt) {

  for (var in vars) {
    new_value <- sid_df[[var]][sid_df$timepoint_w_ == time]
    old_value <- sid_df[[var]][sid_df$timepoint_w_ == rt]
    if (typeof(sid_df[[var]]) != "character") {
      if (absolute_values == "yes") {
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
    if (typeof(sid_df[[var]]) == "character"
    || var == "timepoint_w_" || var == sample_id) {
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

  total_timepoints <- length(unique(df[[input_timepoint]]))

  # use total_timepoints <- 1 to test if not longitudinal
  # If cross sectional, don't make deltas

  # This could probably be simplified
  rf_file_out <- paste0(output_folder, "SELECTED-FEATURES-",
  reference_time, "/", reference_time, ".pdf")

  # TODO this should be reworded, but need to check all downstream uses
  if ((analyze_first != "yes" && reference_time == "first") ||(analyze_previous != "yes" && reference_time == "previous") ||(analyze_pairwise != "yes" && reference_time == "pairwise")) {
    writeLines("Failed Analysis", con = out_file)
    writeLines("Failed Analysis", con = rf_file_out)
    return()
  }

  if (total_timepoints == 1) {
    if (reference_time == "first") {
    writeLines("Failed Analysis", con = out_file)
    writeLines("Failed Analysis", con = rf_file_out)
    return()
    }
  }

  if (total_timepoints <= 2) {
    if (reference_time == "previous" || reference_time == "pairwise") {
      writeLines("Failed Analysis", con = out_file)
      writeLines("Failed Analysis", con = rf_file_out)
      return()
    }
  }

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
  df[[sample_id]] <- as.character(df[[sample_id]])
  # get all variables/features, except ones used for script
  vars <- colnames(df %>% select(!study_id_w_))
  times <- unique(df$timepoint_w_)

  delta_df <- data.frame(sid_delta = character())
  # TODO ref_times are zero in 'previous' improve this.
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
    dplyr::arrange(study_id_w_, timepoint_w_) %>%
    dplyr::select(!sid_delta)

  colnames(delta_df)[colnames(delta_df) == "timepoint_w_"] <- input_timepoint
  colnames(delta_df)[colnames(delta_df) == "study_id_w_"] <- input_study_id

  if (include_reference_values == "no") {
    delta_df <- delta_df %>%
      dplyr::select(-dplyr::ends_with("_reference"))
  }

  delta_df <- eval(parse(text = snakemake@config[["delta_df_mod"]]))

  write.table(delta_df, out_file, row.names = FALSE, sep = "\t")
}


main()