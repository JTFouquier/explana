# Title     : TODO
# Objective : TODO
# Created by: jenniferfouquier
# Created on: 1/12/22

rm(list = ls())
library("tidyverse")
library("faux")

out_file <- "data/simulation/simulation.txt"


my_sample <- function(x, plus_minus) {
  # get a random number plus or minus input value
  my_min <- x - plus_minus
  my_max <- x + plus_minus
  x <- sample(my_min:my_max, 1)[1]
  return(x)
}


list_piece <- function(x, my_factor) {
  # get every "my_factor" or the other items in list
  l1 <- x[seq(1, length(x), my_factor)]
  l2 <- x[!(x %in% l1)]
  return(list(l1, l2))
}


base_dataframe <- function(df, subjects_list, subjects, cohort) {

  df <- df %>%
    filter(study_id %in% subjects_list) %>%
    mutate(cohort = cohort)

  subjects_half1 <- subjects[seq(1, length(subjects), 2)]
  subjects_half2 <- subjects[!(subjects %in% subjects_half1)]


  df <- df %>%
    mutate(therapy = case_when(
      study_id %in% subjects_half1 ~ "therapy_1",
      study_id %in% subjects_half2 ~ "therapy_2"
    ))

  return(df)
}

get_random_intercepts <- function(df, study_id_modifier) {

  total_subjects <- length(unique(df$study_id))

  # random intercepts
  subject_sd <- 3 # SD for the subjects' random intercept

  sub <- tibble(
    study_id = paste0(1:total_subjects, study_id_modifier),
    subject_intercept_effect = rnorm(total_subjects,
    0, subject_sd), # random intercept
  )

  df <- full_join(df, sub)

  return(df)
}


# random slopes
# subject_version_sd <- 20
# subject_i_version_cor <- -0.2
# total_subjects <- 10

# sub <- faux::rnorm_multi(
#   n = total_subjects,
#   vars = 2,
#   r = subject_i_version_cor,
#   mu = 0, # means of random intercepts and slopes are always 0
#   sd = c(subject_sd, subject_version_sd),
#   varnames = c("sub_i", "sub_version_slope")
# ) %>%
#   mutate(
#     study_id = 1:total_subjects,
#     therapy = rep(c("therapy_1", "therapy_2"),
#     each = total_subjects / 2) # between-subjects factor
#   )

# cat_placeholder <- c("", "", "", "", "")
# num_placeholder <- c(0, 0, 0, 0, 0)

c_zeros <- rep(0, 5)
c_char <- rep("", 5)

# Some are placeholders here, but modified later
simulated_df_for_workflow <- function(
  happiness_original = c_zeros,
  cohort = c_char, cohort_effect = c_zeros,
  timepoint = c(1, 2, 3, 4, 5), study_id = c_char,
  microbe_a = c_zeros, microbe_b = c_zeros, microbe_c = c_zeros,
  microbe_a_effect = c_zeros,
  microbe_b_effect = c_zeros,
  microbe_c_effect = c_zeros,
  bonus = rep("no", 5), bonus_effect = c_zeros,
  therapy = c_char, therapy_effect = c_zeros,
  medicine = c_char, medicine_effect = c_zeros,
  income = c_char, income_effect = c_zeros,
  happiness = c_zeros, total_subjects, study_id_modifier
  ) {

  # Create data frame
  study_id <- paste0(1, study_id_modifier)
  study_id <- c(study_id, study_id, study_id, study_id, study_id)
  df <- data.frame(happiness_original, timepoint,
                   cohort, cohort_effect, study_id,
                   microbe_a, microbe_b, microbe_c,
                   microbe_a_effect, microbe_b_effect,
                   microbe_c_effect, bonus, bonus_effect,
                   therapy, therapy_effect, medicine,
                   medicine_effect, income,
                   income_effect, happiness)

  num_columns <- length(colnames(df))

  new_subject_df <- df
  for (i in 2:total_subjects){
      new_subject_df$study_id <- paste0(i, study_id_modifier)
      df <- rbind(df, new_subject_df[rep(1:num_columns, 1)])
  }

  df <- get_random_intercepts(df, study_id_modifier = study_id_modifier)

  subjects <- unique(df$study_id)

  df$sample_id <- paste0(df$study_id, ".", df$timepoint)
  cohort_a_subjects <- subjects[1:(length(subjects) / 2)]
  cohort_b_subjects <- subjects[!subjects %in% cohort_a_subjects]

  df_cohort_a <- base_dataframe(df, cohort_a_subjects, subjects, "A")
  df_cohort_b <- base_dataframe(df, cohort_b_subjects, subjects, "B")
  df <- rbind(df_cohort_a, df_cohort_b)

  # get half of subjects to add variety
  r <- list_piece(subjects, 2)
  subjects_half1 <- r[[1]]
  subjects_half2 <- r[[2]]

  df <- df %>%
    mutate(income =
    case_when(study_id %in% subjects_half1 ~ "high",
              study_id %in% subjects_half2 ~ "low"))

  # Add some no effects to both therapys only if their lifestyle is unhealthy
  # Demonstrates more categorical vars
  df <- df %>%
    group_by(study_id) %>%
    mutate(bonus = case_when(

        # Time 1, 2, 3
        income == "low" & timepoint == 1 ~ "yes",
        income == "low" & timepoint == 2 ~ "yes",
        income == "low" & timepoint == 3 ~ "yes",

        # Time 4
        income == "low" & timepoint == 4
        & cohort == "A" ~ "no",
        income == "low" & timepoint == 4
        & cohort == "B" ~ "yes",

        # Time 5
        income == "low" & timepoint == 5
        & cohort == "B" ~ "no",
        income == "low" & timepoint == 5
        & cohort == "A" ~ "yes",

        income == "high" ~ "yes")) %>%

    mutate(microbe_a = case_when(
      therapy == "therapy_1" & timepoint == 1 ~ my_sample(0.2, 0.03),
      therapy == "therapy_1" & timepoint == 2 ~ my_sample(0.2, 0.03),
      therapy == "therapy_1" & timepoint == 3 ~ my_sample(0.3, 0.03),
      therapy == "therapy_1" & timepoint == 4 ~ my_sample(0.3, 0.03),
      therapy == "therapy_1" & timepoint == 5 ~ my_sample(0.4, 0.03),
      therapy == "therapy_2" ~ 0.2)) %>%

    mutate(microbe_b = case_when(
      (therapy == "therapy_2" & timepoint == 1) ~ 0.1,
      (therapy == "therapy_2" & timepoint == 2) ~ 0.4,
      (therapy == "therapy_2" & timepoint == 3) ~ 0.4,
      (therapy == "therapy_2" & timepoint == 4) ~ 0.4,
      (therapy == "therapy_2" & timepoint == 5) ~ 0.4,
      therapy == "therapy_1" ~ 0.1)) %>%
    ungroup()

  r <- list_piece(subjects, 3)
  med_unaffected <- r[[1]]
  med_affected <- r[[2]]

  r <- list_piece(med_affected, 2)
  med_affected_1 <- r[[1]]
  med_affected_2 <- r[[2]]

  r <- list_piece(med_unaffected, 3)
  med_unaffected_1 <- r[[1]]
  med_unaffected_2 <- r[[2]]


  # all pills have an equal effect on happiness; interaction red_green
  df <- df %>%
  mutate(medicine = case_when(
      study_id %in% med_affected_1 & timepoint %in% 1:2 ~ "GreenPill",
      study_id %in% med_affected_1 & timepoint %in% 3:5 ~ "RedPill",

      study_id %in% med_affected_2 & timepoint %in% 1:2 ~ "RedPill",
      study_id %in% med_affected_2 & timepoint %in% 3:5 ~ "GreenPill",

      # TODO if keeping same pill, change time
      study_id %in% med_unaffected_1 & timepoint %in% 1:2 ~ "RedPill",
      study_id %in% med_unaffected_1 & timepoint %in% 3:5 ~ "RedPill",

      study_id %in% med_unaffected_2 & timepoint %in% 1:2 ~ "GreenPill",
      study_id %in% med_unaffected_2 & timepoint %in% 3:5 ~ "GreenPill",
    ))

  add_effects <- function(df, subjects) {

    for (i in subjects) {
      t1_to_t2 <- subset(df, study_id == i & timepoint %in% 1:2)
      t3_t5 <- subset(df, study_id == i & timepoint %in% 3:5)
      t4 <- subset(df, study_id == i & timepoint == 4)
      t5 <- subset(df, study_id == i & timepoint == 5)

      if ("RedPill" %in% t1_to_t2$medicine
      && "GreenPill" %in% t3_t5$medicine) {
        df$medicine_effect[df$timepoint %in% 3:5
        & df$study_id == i] <- my_sample(-8, 1)
      }

      df$bonus_effect[df$study_id == i & df$bonus == "yes"] <- my_sample(3, 1)

      if ("no" %in% t4$bonus && "yes" %in% t5$bonus) {
        df$bonus_effect[df$timepoint == 5
        & df$study_id == i] <- my_sample(7, 1)
      }

      if ("yes" %in% t4$bonus && "no" %in% t5$bonus) {
        df$bonus_effect[df$timepoint == 5
        & df$study_id == i] <- my_sample(-4, 1)
      }

      # if ("yes" %in% t4$bonus && "yes" %in% t5$bonus) {
      #   df$bonus_effect[df$timepoint == 5
      #   & df$study_id == i] <- my_sample(4, 1)
      # }
    }

  df <- df %>%
    group_by(study_id) %>%
    rowwise() %>%
    mutate(cohort_effect = case_when(
           cohort == "A" ~ my_sample(-6, 1),
           cohort == "B" ~ my_sample(6, 1)),

           microbe_b_effect = microbe_b * 4,

           income_effect = case_when(
           income == "high" ~ my_sample(3, 1),
           income == "low" ~ my_sample(-3, 1)),

           microbe_a_effect = microbe_a * 2) %>%
    ungroup()

  df_i1 <- df %>%
    filter(therapy == "therapy_1") %>%
    group_by(study_id) %>%
    rowwise() %>%
    mutate(therapy_effect = case_when(
           timepoint == 1 ~ my_sample(0, 0),
           timepoint == 2 ~ my_sample(5, 1),
           timepoint == 3 ~ my_sample(10, 1),
           timepoint == 4 ~ my_sample(20, 1),
           timepoint == 5 ~ my_sample(40, 1))) %>%
    ungroup()

  df_i2 <- df %>%
    filter(therapy == "therapy_2") %>%
    group_by(study_id) %>%
    rowwise() %>%
    mutate(therapy_effect = case_when(
           timepoint == 1 ~ my_sample(0, 0),
           timepoint == 2 ~ my_sample(4, 1),
           timepoint == 3 ~ my_sample(8, 1),
           timepoint == 4 ~ my_sample(12, 1),
           timepoint == 5 ~ my_sample(16, 1))) %>%
    ungroup()

  # merge both dataframes back together and sort
  df <- rbind(df_i1, df_i2)
  df <- df %>% arrange(study_id)

  df <- df %>%
    mutate(happiness = happiness_original +
    rowSums(.[grep("_effect", names(.))]))

  write.table(df, "data/simulation/simulation-effects.txt",
  row.names = FALSE, sep = "\t")
  df <- df %>%
    select(-ends_with("_effect")) %>%
    select(!(happiness_original))
  return(df)
  }

  # remove effect columns and original happiness
  df <- add_effects(df, subjects)
  return(df)
}

df_a <- simulated_df_for_workflow(happiness_original = rep(30, 5),
cohort = rep("A", 5),
total_subjects = 10, study_id_modifier = "a")

df_b <- simulated_df_for_workflow(happiness_original = rep(30, 5),
cohort = rep("B", 5),
total_subjects = 10, study_id_modifier = "b")
print(df_a)
print(df_b)
df <- rbind(df_a, df_b)
print(df)

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

  file_name <- paste0(file_prefix, num_features, "-random-vars.txt")
  write.table(df, file_name, row.names = FALSE, sep = "\t")

  file_name <- paste0(file_prefix, num_features,
  "-random-vars-shuffled.txt")
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

addFeatures("data/simulation/simulation.txt", "happiness",
            "data/simulation/simulation-", 10, 100, 1000)
