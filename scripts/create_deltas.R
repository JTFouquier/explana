
library(dplyr)
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(usedist)
source("scripts/viz-datatable.R")

in_file = snakemake@input[["in_file"]]
out_file = snakemake@output[["out_file"]]
output_dir = snakemake@config[["out"]]

reference_time = snakemake@params[["reference_time"]]
absolute_values = snakemake@params[["absolute_values"]]
build_visualizer = snakemake@params[["build_visualizer"]]

# TODO allow a list of dms
dm_file = snakemake@params[["distance_matrix"]]
dm = read.table(dm_file, fill=T)

# TODO add option for looking at only subset of times
#times_evaluated = snakemake@params[["times_evaluated"]]

df = read.csv(in_file, sep = "\t", check.names = FALSE)

df$Timepoint = factor(df$Timepoint, ordered = TRUE)
df$StudyID.Timepoint = NULL
# get all variables/features, except ones used for script
#vars <- colnames(df %>% select(!c(Timepoint, StudyID)))
vars <- colnames(df %>% select(!StudyID))

times = unique(df$Timepoint)

delta.df = data.frame(StudyID.Timepoint = character())
# TODO ref.times are zero in 'previous' improve this.
# TODO make function for ref.times
# TODO add StudyID.Timepoint AKA SampleID

# TODO add ability to use .QZA or TSV
#unweighted.unifrac.file = "data/unweighted_unifrac_distance_matrix.qza"
#dm <- read_qza(unweighted.unifrac.file)
#dm = dm$data

# TODO allow multiple distance matrix files to be added via function

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
      studyid.df <- df %>% filter((Timepoint == rt | Timepoint == time) & StudyID == studyid)
      # skip if no samples for subject at time points of interest
      # also skip if reference time is larger than time (would be duplicates)
      if (nrow(studyid.df) < 2 || rt > time){ next }
      dm_value = dist_subset(dm, c(paste0(studyid, ".", rt),
                                 paste0(studyid, ".", time)))[1]
      new.comparison.id = paste0(studyid, "_", rt, "_", time)
      df.new = data.frame(StudyID.Timepoint = new.comparison.id,
                          Timepoint = paste0(rt, "_", time), # TODO needed?
                          StudyID = studyid,
                          dist_matrix = dm_value)
      delta.df = dplyr::bind_rows(df.new, delta.df)

      for (var in vars) {
        new.value = studyid.df[[var]][studyid.df$Timepoint == time]
        old.value = studyid.df[[var]][studyid.df$Timepoint == rt]
        if (typeof(studyid.df[[var]]) != "character"){
          if (reference_time == "pairwise" && absolute_values == "yes"){
            # TODO convert to abs values at end? reduce computational time
            delta.df[[var]][delta.df$StudyID.Timepoint == new.comparison.id] = abs(new.value - old.value)
            delta.df[[paste0(var,"_reference")]][delta.df$StudyID.Timepoint == new.comparison.id] = old.value
          } else {
            delta.df[[var]][delta.df$StudyID.Timepoint == new.comparison.id] = new.value - old.value
            delta.df[[paste0(var,"_reference")]][delta.df$StudyID.Timepoint == new.comparison.id] = old.value
          }
        }
        if (typeof(studyid.df[[var]]) == "character" || var == "Timepoint"){
          delta.df[[var]][delta.df$StudyID.Timepoint == new.comparison.id] = paste0(old.value, "_", new.value)
          delta.df[[paste0(var,"_reference")]][delta.df$StudyID.Timepoint == new.comparison.id] = paste0(old.value)
        }
        else{ next }
      }
    }
  }
}

delta.df <- delta.df %>% arrange(StudyID.Timepoint, Timepoint)
write.table(delta.df, out_file, row.names = FALSE, sep = "\t")

# TODO optional create visualizer
# TODO change this to either df or load the dataframe
if (build_visualizer == TRUE){
  input_file_name = paste0(output_dir, "04-SELECTED-FEATURES-", reference_time,
                           "/deltas-", reference_time , ".txt")
  output_file_name = paste0(output_dir, "04-SELECTED-FEATURES-",
                            reference_time, "/vizualizer-", reference_time,
                            ".html")
  build_datatable(input_file_name, output_file_name)
}