
library(dplyr)
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(usedist)
source("scripts/viz-datatable.R")

# access data from Snakemake's Snakefile rules and user's config.yaml
in_file = snakemake@input[["in_file"]]
out_file = snakemake@output[["out_file"]]
output_folder = snakemake@config[["out"]]
sample_id = snakemake@config[["sample_id"]]

reference_time = snakemake@params[["reference_time"]]
absolute_values = snakemake@params[["absolute_values"]]
build_datatable = snakemake@params[["build_datatable"]]

# get name of distance matrix in params for rule, then
distance_matrices = eval(parse(text = snakemake@params[["distance_matrices"]]))
distance_matrices <- lapply(distance_matrices, read.table, fill=T)

if (length(distance_matrices) > 0){
  dm_flag = TRUE
}

df = read.csv(in_file, sep = "\t", check.names = FALSE)

df$Timepoint = factor(df$Timepoint, ordered = TRUE)
df$StudyID.delta = NULL
# get all variables/features, except ones used for script
vars <- colnames(df %>% select(!StudyID))
times = unique(df$Timepoint)

delta.df = data.frame(StudyID.delta = character())
# TODO ref.times are zero in 'previous' improve this.
# TODO make function for ref.times

# TODO add ability to use .QZA or TSV
#unweighted.unifrac.file = "data/unweighted_unifrac_distance_matrix.qza"
#dm <- read_qza(unweighted.unifrac.file)
#dm = dm$data


for (studyid in unique(df$StudyID)){
  for (time in times) {
    time = as.numeric(time)
    # get different reference times to create deltas
    if (reference_time == "first"){
      ref.times = 1
    } else if (reference_time == "previous"){
      ref.times = time - 1
    }
    else if (reference_time == "pairwise"){
      ref.times = as.numeric(levels(df$Timepoint)) # TODO times or Timepoint
      ref.times = ref.times[ref.times != time]
    }
    for (rt in ref.times){
      studyid.df <- df %>%
        filter((Timepoint == rt | Timepoint == time) & StudyID == studyid)
      # skip if no samples for subject at time points of interest
      # also skip if reference time is larger than time (would be duplicates)
      if (nrow(studyid.df) < 2 || rt > time){ next }
      ref_sample = studyid.df %>%
        filter(StudyID == {{studyid}} & Timepoint == {{rt}})
      current_sample = studyid.df %>%
        filter(StudyID == {{studyid}} & Timepoint == {{time}})
      comparison = paste0(studyid, "_", rt, "_", time)
      # Create new dataframe for each comparison
      df.new.comparison = data.frame(StudyID.delta = comparison,
                                     Timepoint = paste0(rt, "_", time),
                                     StudyID = studyid)
      # add all distances (comparisons between sample timepoints); NA if none
      # from list of distance matrix files
      if (dm_flag){
        for (i in names(distance_matrices)){
          result = tryCatch({
            result = dist_subset(distance_matrices[[i]], c(ref_sample$SampleID,
                                 current_sample$SampleID))[1]
          }, error = function(err){
            result = NA
          }, finally = function(){
            return(result)
          })
        df.new.comparison = df.new.comparison %>%
          mutate({{i}} := result)
        }
      }
      delta.df = dplyr::bind_rows(df.new.comparison, delta.df)
      for (var in vars) {
        new.value = studyid.df[[var]][studyid.df$Timepoint == time]
        old.value = studyid.df[[var]][studyid.df$Timepoint == rt]
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
        if (typeof(studyid.df[[var]]) == "character" || var == "Timepoint"){
          delta.df[[var]][delta.df$StudyID.delta == comparison] =
            paste0(old.value, "_", new.value)
          delta.df[[paste0(var,"_reference")]][delta.df$StudyID.delta == comparison] = paste0(old.value)
        }
        else{ next }
      }
    }
  }
}

delta.df <- delta.df %>%
  arrange(StudyID, Timepoint) %>%
  select(!StudyID.delta)

write.table(delta.df, out_file, row.names = FALSE, sep = "\t")

# TODO optional create visualizer
# TODO change this to either df or load the dataframe
if (build_datatable %in% c("TRUE", "True")){
  input_file_name = paste0(output_folder, "04-SELECTED-FEATURES-", reference_time,
                           "/deltas-", reference_time , ".txt")
  output_file_name = paste0(output_folder, "04-SELECTED-FEATURES-",
                            reference_time, "/vizualizer-", reference_time,
                            ".html")
  build_datatable_viz(input_file_name, output_file_name)
}