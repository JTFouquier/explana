# Title     : TODO
# Objective : TODO
# Created by: jenniferfouquier
# Created on: 1/12/22

rm(list = ls())

library("tidyverse")

total_subjects <- 80
out_file <- "data/simulated-scripted-test-04032023.txt"

  # add random features
#low <- 50
#med <- 100
#high <- 200
### 1_2 (1) 1_3 (2)
## 1_2 (1) 2_3 (1)


# Some are placeholders here, but modified later
inflammation_original <- c(25, 25, 25)
cohort <- c("", "", "")
cohort_effect <- c(0, 0, 0)
timepoint <- c(1, 2, 3)
study_id <- c(1, 1, 1)
microbe_a <- c(0, 0, 0)
microbe_b <- c(0.10, 0.15, 0.15)
adhered <- c("Yes", "Yes", "Yes")
adhered_effect <- c(0, 0, 0)
diet <- c("", "", "")
diet_effect <- c(0, 0, 0)
medicine <- c("", "", "")
medicine_effect <- c(0, 0, 0)
healthy_lifestyle <- c("", "", "")
healthy_lifestyle_effect <- c(0, 0, 0)
inflammation <- c(0, 0, 0)

# Create data frame
df <- data.frame(inflammation_original, timepoint, cohort, cohort_effect,
                 study_id, microbe_a, microbe_b, adhered, adhered_effect,
                 diet, diet_effect, medicine, medicine_effect,
                 healthy_lifestyle, healthy_lifestyle_effect, inflammation)

num_columns <- length(colnames(df))

new_subject_df <- df
for (i in 2:total_subjects){
    new_subject_df$study_id <- i
    df <- rbind(df, new_subject_df[rep(1:num_columns, 1)])
}

subjects <- unique(df$study_id)

df$sample_id <- paste0(df$study_id, ".", df$timepoint)
overweight_subjects <- subjects[1:(length(subjects) / 2)]
healthy_subjects <- subjects[!subjects %in% overweight_subjects]


base_dataframe <- function(df, subjects_list, cohort) {

  cohort_df <- df %>%
    filter(study_id %in% subjects_list) %>%
    mutate(cohort = cohort)

  h1 <- subjects_list[1:(length(subjects_list) / 2)]
  h2 <- subjects_list[!subjects_list %in% h1]

  cohort_df <- cohort_df %>%
    mutate(diet = case_when(
      study_id %in% h1 ~ "diet_1",
      study_id %in% h2 ~ "diet_2"
    ))

  return(cohort_df)
}
 
df_overweight <- base_dataframe(df, overweight_subjects, "overweight")
df_healthy <- base_dataframe(df, healthy_subjects, "healthy")
df <- rbind(df_overweight, df_healthy)


# healthylifestyle odd & unhealthylifestyle even
odd_subjects <- subjects[seq(1, length(subjects), 2)]
even_subjects <- subjects[!(subjects %in% odd_subjects)]

df <- df %>%
  mutate(healthy_lifestyle = case_when(study_id %in% odd_subjects ~ "Yes",
                                       study_id %in% even_subjects ~ "No"))

# Add some no effects to both diets only if their lifestyle is unhealthy
# Demonstrates more categorical vars
df <- df %>%
  mutate(adhered = case_when(
      healthy_lifestyle == "No" & timepoint == 2 & diet == "diet_1" ~ "No",
      healthy_lifestyle == "No" & timepoint == 2 & diet == "diet_2" ~ "Yes",
      healthy_lifestyle == "No" & timepoint == 3 & diet == "diet_2" ~ "No",
      healthy_lifestyle == "No" & timepoint == 3 & diet == "diet_1" ~ "Yes",
      healthy_lifestyle == "No" & timepoint == 1 ~ "Yes",
      healthy_lifestyle == "Yes" ~ "Yes")) %>%
  mutate(microbe_a = case_when(
    diet == "diet_1" & timepoint == 1 ~ 0.1,
    diet == "diet_1" & timepoint == 2 ~ 0.2,
    diet == "diet_1" & timepoint == 3 ~ 0.4,
    diet == "diet_2" ~ 0.1))

# every x subject in first list then remaining in other list
med_affected <- subjects[seq(1, length(subjects), 4)]
med_unaffected <- subjects[!(subjects %in% med_affected)]

med_affected_few <- med_affected[seq(1, length(med_affected), 4)]
med_affected_most <- med_affected[!(med_affected %in% med_affected_few)]

df <- df %>%
  mutate(medicine = case_when(
      #study_id %in% med_unaffected & (Timepoint == 1 | Timepoint == 2) ~ "G",
      #study_id %in% med_unaffected & Timepoint == 3 ~ "P",
      study_id %in% med_affected_few & timepoint == 1 ~ "G", ### F A G
      study_id %in% med_affected_few & timepoint == 2 ~ "A",
      study_id %in% med_affected_few & timepoint == 3 ~ "F",
      #study_id %in% med_affected_most & Timepoint == 1 ~ "P",
      study_id %in% med_affected_most & timepoint == 2 ~ "F", ###
      study_id %in% med_affected_most & timepoint == 3 ~ "G" ###
))


meds <- c("G", "A", "F")
r <- length(df$medicine[df$study_id %in% med_unaffected])
df$medicine[df$study_id %in% med_unaffected] <- sample(meds, r, replace = TRUE)
r <- length(df$medicine[df$study_id %in% med_affected_most & timepoint == 1])
df$medicine[df$study_id %in% med_affected_most & timepoint == 1] <- sample(meds, r, replace = TRUE)

# ADD EFFECTS
for (i in subjects){
  df_subject_1_2 <- subset(df, study_id == i & !timepoint == 3)
  df_subject_3 <- subset(df, study_id == i & timepoint == 3)
  if ("F" %in% df_subject_1_2$medicine && "G" %in% df_subject_3$medicine) {
    df$medicine_effect[df$timepoint == 3 & df$study_id == i] <- 10
  }
}

df$cohort_effect[df$cohort == "Overweight"] <- 25
df$adhered_effect[df$adhered == "No"] <- 4
df$healthy_lifestyle_effect[df$healthy_lifestyle == "Yes"] <- -3

df$diet_effect[df$timepoint == 2 & df$diet == "diet_1"] <- -20
df$diet_effect[df$timepoint == 3 & df$diet == "diet_1"] <- -30

df$diet_effect[df$timepoint == 2 & df$diet == "diet_2"] <- -5
df$diet_effect[df$timepoint == 3 & df$diet == "diet_2"] <- -6


df <- df %>%
  mutate(inflammation = inflammation_original +
  rowSums(.[grep("_effect", names(.))]))

write.table(df, "data/simulated-scripted-test-effects-04032023.txt",
row.names = FALSE, sep = "\t")
# remove effect columns and original inflammation
df <- df %>%
  select(-ends_with("_effect")) %>%
  select(!(inflammation_original))

## Calculate output
write.table(df, out_file, row.names = FALSE, sep = "\t")





writeFiles <- function(file_prefix, df, num_features, num_features_prev,
                       response_var, sampled) {

  for (i in num_features_prev + 1:num_features){
    feature_name <- paste0('V', i)
    df[[feature_name]] <- sample(1000, size = nrow(df), replace = TRUE)
  }

  df_shuffled_response <- df
  df_shuffled_response[[response_var]] <- sampled

  file_name <- paste0(file_prefix, num_features, "-random-vars-04032023.txt")
  write.table(df, file_name, row.names = FALSE, sep = "\t")

  file_name <- paste0(file_prefix, num_features, "-random-vars-shuffled-04032023.txt")
  write.table(df_shuffled_response, file_name, row.names = FALSE, sep = "\t")

  return(df)
}


# add shuffle response option
addFeatures <- function(df, response_var, file_prefix, low, med, hi) {

  df <- read.csv(df, sep = "\t")

  # shuffle response
  shuffled_response <- sample(df[[response_var]])

  # TODO this is repetitive (make list starting with 0)
  df <- writeFiles(file_prefix, df, low, 0, response_var,
                   shuffled_response)
  df <- writeFiles(file_prefix, df, med, low, response_var, shuffled_response)
  df <- writeFiles(file_prefix, df, hi, med, response_var, shuffled_response)

}

addFeatures("data/simulated-scripted-test-04032023.txt", "inflammation",
            "data/simulated-scripted-test-04032023.txt", 10, 100, 1000)
