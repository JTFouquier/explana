# Title     : TODO
# Objective : TODO
# Created by: jenniferfouquier
# Created on: 1/12/22

rm(list=ls())

library("tidyverse")

total_subjects = 24
out_file = "data/simulated-scripted-TEST.txt"

  # add random features
low = 50
med = 100
high = 200


# Some are placeholders here, but modified later
Inflammation_original = c(25, 25, 25)
Cohort = c("", "", "")
Cohort_effect = c(0, 0, 0)
Timepoint = c(1, 2, 3)
StudyID = c(1, 1, 1)
MicrobeA = c(0, 0, 0)
MicrobeB = c(0.10, 0.15, 0.15)
Adhered = c("Yes", "Yes", "Yes")
Adhered_effect = c(0, 0, 0)
Diet = c("", "", "")
Diet_effect = c(0, 0, 0)
Medicine = c("", "", "")
Medicine_effect = c(0, 0, 0)
HealthyLifestyle = c("", "", "")
HealthyLifestyle_effect = c(0, 0, 0)
Inflammation = c(0, 0, 0)

# Create data frame
df <- data.frame(Inflammation_original, Timepoint, Cohort, Cohort_effect,
                 StudyID, MicrobeA, MicrobeB, Adhered, Adhered_effect,
                 Diet, Diet_effect, Medicine, Medicine_effect,
                 HealthyLifestyle, HealthyLifestyle_effect, Inflammation)

num_columns = length(colnames(df))

new_subject_df = df
for (i in 2:total_subjects){
    new_subject_df$StudyID = i
    df = rbind(df, new_subject_df[rep(1:num_columns, 1)])
}

subjects = unique(df$StudyID)
overweight_subjects = subjects[1:(length(subjects)/2)]
healthy_subjects = subjects[!subjects %in% overweight_subjects]


base_dataframe = function(df, subjects_list, cohort){

  cohort_df = df %>%
    filter(StudyID %in% subjects_list) %>%
    mutate(Cohort = cohort)

  h1 = subjects_list[1:(length(subjects_list)/2)]
  h2 = subjects_list[!subjects_list %in% h1]

  cohort_df = cohort_df %>%
    mutate(Diet = case_when(
      StudyID %in% h1 ~ "Keto",
      StudyID %in% h2 ~ "Western"
    ))

  return(cohort_df)
}

df_overweight = base_dataframe(df, overweight_subjects, "Overweight")
df_healthy = base_dataframe(df, healthy_subjects, "Healthy")
df = rbind(df_overweight, df_healthy)


# healthylifestyle odd & unhealthylifestyle even
odd_subjects = subjects[seq(1,length(subjects),2)]
even_subjects = subjects[!(subjects %in% odd_subjects)]

df = df %>%
  mutate(HealthyLifestyle = case_when(StudyID %in% odd_subjects ~ "Yes",
                                      StudyID %in% even_subjects ~ "No"))

df = df %>%
  mutate(Adhered = case_when(
      HealthyLifestyle == "No" & Timepoint == 2 & Diet == "Keto" ~ "No",
      HealthyLifestyle == "No" & Timepoint == 2 & Diet == "Western" ~ "Yes",
      HealthyLifestyle == "No" & Timepoint == 3 & Diet == "Western" ~ "No",
      HealthyLifestyle == "No" & Timepoint == 3 & Diet == "Keto" ~ "Yes",
      HealthyLifestyle == "No" & Timepoint == 1 ~ "Yes",
      HealthyLifestyle == "Yes" ~ "Yes")) %>%
  mutate(MicrobeA = case_when(
    Diet == "Keto" & Timepoint == 1 ~ 0.1,
    Diet == "Keto" & Timepoint == 2 ~ 0.2,
    Diet == "Keto" & Timepoint == 3 ~ 0.4,
    Diet == "Western" ~ 0.0))

s1 = subjects[seq(1,length(subjects),5)]
s2 = subjects[!(subjects %in% s1)]

s1_few = s1[seq(1,length(s1),5)]
s1_most = s1[!(s1 %in% s1_few)]

df = df %>%
  mutate(Medicine = case_when(
      StudyID %in% s2 & (Timepoint == 1 | Timepoint == 2) ~ "G",
      StudyID %in% s2 & Timepoint == 3 ~ "F",
      StudyID %in% s1_few & Timepoint == 1 ~ "F",
      StudyID %in% s1_few & Timepoint == 2 ~ "P",
      StudyID %in% s1_few & Timepoint == 3 ~ "G",
      StudyID %in% s1_most & Timepoint == 1 ~ "P",
      StudyID %in% s1_most & Timepoint == 2 ~ "F",
      StudyID %in% s1_most & Timepoint == 3 ~ "G"
))

# ADD EFFECTS
for (i in subjects){
  df_subject_1_2 = subset(df, StudyID == i & !Timepoint == 3)
  df_subject_3 = subset(df, StudyID == i & Timepoint == 3)
  if ("F" %in% df_subject_1_2$Medicine & "G" %in% df_subject_3$Medicine){
    df$Medicine_effect[df$Timepoint==3 & df$StudyID==i] = 7
  }
}

df$Cohort_effect[df$Cohort == "Overweight"] = 25
df$Adhered_effect[df$Adhered == "No"] = 4
df$HealthyLifestyle_effect[df$HealthyLifestyle == "Yes"] = -3

df$Diet_effect[df$Timepoint == 2 & df$Diet == "Keto"] = -20
df$Diet_effect[df$Timepoint == 3 & df$Diet == "Keto"] = -30

df$Diet_effect[df$Timepoint == 2 & df$Diet == "Western"] = -5
df$Diet_effect[df$Timepoint == 3 & df$Diet == "Western"] = -6


df = df %>%
  mutate(Inflammation = Inflammation_original + rowSums(.[grep("_effect", names(.))]))

# remove effect columns and original inflammation
df = df %>%
  select(-ends_with("_effect")) %>%
  select(!(Inflammation_original))

## Calculate output
write.table(df, out_file, row.names = FALSE, sep = "\t")


writeFiles = function(file_prefix, df, num_features, num_features_prev,
                      response_var, sampled){

  for (i in num_features_prev + 1:num_features){
    feature_name = paste0('V', i)
    df[[feature_name]] <- sample(1000, size = nrow(df), replace = TRUE)
  }

  df_shuffled_response = df
  df_shuffled_response[[response_var]] = sampled

  file_name = paste0(file_prefix, num_features, "-random-vars.txt")
  write.table(df, file_name, row.names = FALSE, sep = "\t")

  file_name = paste0(file_prefix, num_features, "-random-vars-shuffled.txt")
  write.table(df_shuffled_response, file_name, row.names = FALSE, sep = "\t")

  return(df)
}


# add shuffle response option
addFeatures <- function(df, response_var, file_prefix, low, med, hi){

  df = read.csv(df, sep = "\t")

  # shuffle response
  shuffled_response = sample(df[[response_var]])

  # TODO this is repetitive
  df = writeFiles(file_prefix, df, low, 0, response_var,
                  shuffled_response)
  df = writeFiles(file_prefix, df, med, low, response_var, shuffled_response)
  df = writeFiles(file_prefix, df, hi, med, response_var, shuffled_response)

}

addFeatures("data/simulated-scripted-TEST.txt", "Inflammation",
            "data/simulated-scripted-TEST-", 50, 100, 200)




