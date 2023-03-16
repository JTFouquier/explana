
#devtools::install_github("jbisanz/qiime2R")
# library(qiime2R)

source("scripts/install.R")
package_list <- c("dplyr", "usedist", "utils")
install_r_packages(package_list = package_list)

library(dplyr)
library(usedist)
library(utils)

source("scripts/viz-datatable.R")

# from config
output_folder = snakemake@config[["out"]]

sample_id = snakemake@config[["sample_id"]]

# change back to these names
input_timepoint = snakemake@config[["timepoint"]]
input_study_id = snakemake@config[["random_effect"]]


# from input, output, and params
in_file = snakemake@input[["in_file"]]
out_file = snakemake@output[["out_file"]]

reference_time = snakemake@params[["reference_time"]]
absolute_values = snakemake@params[["absolute_values"]]
build_datatable = snakemake@params[["build_datatable"]]

# get name of distance matrix in params for rule, then
distance_matrices = eval(parse(text = snakemake@params[["distance_matrices"]]))
distance_matrices <- lapply(distance_matrices, read.table, fill=T)


include_distance_matrices = function(distance_matrices, ref_sample,
                                     current_sample, df.new.comparison){
  # add all distances (comparisons between sample timepoints); NA if none
    for (i in names(distance_matrices)){
      result = tryCatch({
        result =
          dist_subset(distance_matrices[[i]],
                      c(ref_sample[[sample_id]],
                        current_sample[[sample_id]]))[1]
      }, error = function(err){
        result = NA
      }, finally = function(){
        return(result)
      })
      df.new.comparison = df.new.comparison %>%
        mutate({{i}} := result)
    }
  return(df.new.comparison) # update the new comparison
}



diffs_for_all_vars_per_subject = function(delta.df, vars, studyid.df,
                                          comparison, time, rt){

  for (var in vars) {
    new.value = studyid.df[[var]][studyid.df$Timepoint_w_ == time]
    old.value = studyid.df[[var]][studyid.df$Timepoint_w_ == rt]
    if (typeof(studyid.df[[var]]) != "character"){
      if (reference_time == "pairwise" && absolute_values == "yes"){
        # TODO convert to abs values at end? reduce computational time
        delta.df[[var]][delta.df$StudyID.delta == comparison] = abs(new.value - old.value)
        delta.df[[paste0(var,"_reference")]][delta.df$StudyID.delta == comparison] = old.value
      } else {
        delta.df[[var]][delta.df$StudyID.delta == comparison] = new.value - old.value
        delta.df[[paste0(var,"_reference")]][delta.df$StudyID.delta == comparison] = old.value
      }
    }
    if (typeof(studyid.df[[var]]) == "character" || var == "Timepoint_w_"){
      delta.df[[var]][delta.df$StudyID.delta == comparison] =
        paste0(old.value, "_", new.value)
      delta.df[[paste0(var,"_reference")]][delta.df$StudyID.delta == comparison] = paste0(old.value)
    }
    else{ next }
  }

  return(delta.df)
}


all_reference_time_comparisons = function(df, studyid, ref.times, time, vars,
                                          delta.df, dm_flag){
  for (rt in ref.times){
    studyid.df <- df %>%
      filter((Timepoint_w_ == rt | Timepoint_w_ == time) & StudyID_w_ == studyid)
    # skip if no samples for subject at time points of interest
    # also skip if reference time is larger than time (would be duplicates)
    if (nrow(studyid.df) < 2 || rt > time){ next }
    ref_sample = studyid.df %>%
      filter(StudyID_w_ == {{studyid}} & Timepoint_w_ == {{rt}})
    current_sample = studyid.df %>%
      filter(StudyID_w_ == {{studyid}} & Timepoint_w_ == {{time}})
    comparison = paste0(studyid, "_", rt, "_", time)

    # Create new dataframe for each comparison
    df.new.comparison = data.frame(StudyID.delta = comparison,
                                   Timepoint_w_ = paste0(rt, "_", time),
                                   StudyID_w_ = studyid)

    # if there are dms, add differences here for each subject
    if (dm_flag){
      df.new.comparison =
        include_distance_matrices(distance_matrices = distance_matrices,
                                  ref_sample = ref_sample,
                                  current_sample = current_sample,
                                  df.new.comparison = df.new.comparison)
    }

    delta.df = dplyr::bind_rows(df.new.comparison, delta.df)

    # calculate differences for all variables for one subject
    delta.df =
      diffs_for_all_vars_per_subject(delta.df = delta.df, vars = vars,
                                     studyid.df = studyid.df,
                                     comparison = comparison,
                                     time = time, rt = rt)
  }
  return(delta.df)
}


main = function(){
  dm_flag = FALSE
  if (length(distance_matrices) > 0){
    dm_flag = TRUE
  }

  df = read.csv(in_file, sep = "\t", check.names = FALSE)

  # if "Timepoint_w_" or "StudyID_w_" exist here, then can't go on
  if ("Timepoint_w_" %in% colnames(df)) { print("Item is present in the List.")}
  if ("StudyID_w_" %in% colnames(df)){ print("Item is present in the List.")}

  colnames(df)[colnames(df) == input_timepoint] <- "Timepoint_w_"
  colnames(df)[colnames(df) == input_study_id] <- "StudyID_w_"

  df$Timepoint_w_ = factor(df$Timepoint_w_, ordered = TRUE)
  df$StudyID.delta = NULL
  # get all variables/features, except ones used for script
  vars <- colnames(df %>% select(!StudyID_w_))
  times = unique(df$Timepoint_w_)

  delta.df = data.frame(StudyID.delta = character())
  # TODO ref.times are zero in 'previous' improve this.
  # TODO make function for ref.times

  # TODO add ability to use .QZA or TSV
  #unweighted.unifrac.file = "data/unweighted_unifrac_distance_matrix.qza"
  #dm <- read_qza(unweighted.unifrac.file)
  #dm = dm$data

  # calculate all reference time comparisons, per var, per person

  

  for (studyid in unique(df$StudyID_w_)){
    for (time in times) {
      time = as.numeric(time)
      # get different reference times to create deltas
      if (reference_time == "first"){
        ref.times = 1
      } else if (reference_time == "previous"){
        ref.times = time - 1
      }
      else if (reference_time == "pairwise"){
        ref.times = as.numeric(levels(df$Timepoint_w_)) # TODO times or Timepoint
        ref.times = ref.times[ref.times != time]
      }

      delta.df = all_reference_time_comparisons(df = df, studyid = studyid,
                                                ref.times = ref.times,
                                                time = time,
                                                vars = vars,
                                                delta.df = delta.df,
                                                dm_flag = dm_flag)
    }
  }

  delta.df <- delta.df %>%
    arrange(StudyID_w_, Timepoint_w_) %>%
    select(!StudyID.delta)

  colnames(delta.df)[colnames(delta.df) == "Timepoint_w_"] <- input_timepoint
  colnames(delta.df)[colnames(delta.df) == "StudyID_w_"] <- input_study_id

  write.table(delta.df, out_file, row.names = FALSE, sep = "\t")

  # TODO optional create visualizer
  # TODO change this to either df or load the dataframe
  if (build_datatable %in% c("TRUE", "True")){
    input_file_name = paste0(output_folder, "04-SELECTED-FEATURES-", reference_time,
                             "/", reference_time , ".txt")
    output_file_name = paste0(output_folder, "04-SELECTED-FEATURES-",
                              reference_time, "/vizualizer-", reference_time,
                              ".html")
    build_datatable_viz(input_file_name, output_file_name)
  }
}


main()