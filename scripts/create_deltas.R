
library(dplyr)

in_file = snakemake@input[["in_file"]]
out_file = snakemake@output[["out_file"]]

reference_time = snakemake@params[["reference_time"]]
absolute_values = snakemake@params[["absolute_values"]]
# TODO add option for looking at only subset of times
#times_evaluated = snakemake@params[["times_evaluated"]]

df = read.csv(in_file, sep = "\t")

df$Timepoint = factor(df$Timepoint, ordered = TRUE)

# get all variables/features, except ones used for script
vars <- colnames(df %>% select(!c(Timepoint, StudyID)))

times = unique(df$Timepoint)

delta.df = data.frame(StudyID.Timepoint = character())

# TODO ref.times are zero in 'previous' improve this.
# TODO make function for ref.times
# TODO add StudyID.Timepoint AKA SampleID

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

      new.comparison.id = paste0(studyid, "_", rt, "_", time)
      df.new = data.frame(StudyID.Timepoint = new.comparison.id,
                          Timepoint = paste0(rt, "_", time),
                          StudyID = studyid)
      delta.df = dplyr::bind_rows(df.new, delta.df)

      for (var in vars) {
        new.value = studyid.df[[var]][studyid.df$Timepoint == time]
        old.value = studyid.df[[var]][studyid.df$Timepoint == rt]
        if (typeof(studyid.df[[var]]) != "character"){

          if (reference_time == "pairwise" && absolute_values == "yes"){
            # TODO convert to abs values at end? reduce computational time
            delta.df[[var]][delta.df$StudyID.Timepoint == new.comparison.id] = abs(new.value - old.value)
          } else {
            delta.df[[var]][delta.df$StudyID.Timepoint == new.comparison.id] = new.value - old.value
          }
        }
        else{
          delta.df[[var]][delta.df$StudyID.Timepoint == new.comparison.id] = paste0(old.value, "_", new.value)
        }
      }
    }
  }
}

delta.df <- delta.df %>% arrange(StudyID.Timepoint, Timepoint)
write.table(delta.df, out_file, row.names = FALSE, sep = "\t")