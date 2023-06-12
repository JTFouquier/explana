# Title     : TODO
# Objective : TODO
# Created by: jenniferfouquier
# Created on: 1/12/22

rm(list = ls())
library("tidyverse")
library("faux")

out_file <- "data/simulation/simulation.txt"


# vars

# two cohorts (A and)
subjects_per_cohort <- 40
happiness_start <- 60

std_dev_random_intercepts <- 7
total_timepoints <- 5

# for initializing dataframe
c_zeros <- rep(0, total_timepoints)
c_char <- rep("", total_timepoints)

green_blue_effect <- -15

cohort_a_effect <- -9
cohort_b_effect <- -4

income_high_effect <- 5
income_low_effect <- -5

relationship_s_effect <- -4
relationship_p_effect <- 4
relationship_s_p_effect <- 6
relationship_p_s_effect <- -10

# Explanation of variables

# microbe_b (numerical)
# microbe b - starts low but levels out after first timepoint
# in all individuals

# microbe_a (numerical)
# increases in therapy_1, stays the same therapy_2

# medication (categorical)
# all pills have an equal effect on happiness; interaction Green_Blue causes negative impact on happiness; blue is not found at timepoint 1 to show that this interaction is unable to be found using the First delta dataset. Change is at end of study. 

# relationship (categorical)
# shows important in first deltas (change is at start)

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
  subject_sd <- std_dev_random_intercepts

  sub <- tibble(
    study_id = paste0(1:total_subjects, study_id_modifier),
    subject_intercept_effect = rnorm(total_subjects,
    0, subject_sd), # random intercept
  )

  df <- full_join(df, sub)

  return(df)
}

# Some are placeholders here, but modified later
simulated_df_for_workflow <- function(
  happiness_original = c_zeros,
  cohort = c_char, cohort_effect = c_zeros,
  timepoint = c(1, 2, 3, 4, 5), study_id = c_char,
  microbe_a = c_zeros, microbe_b = c_zeros,
  microbe_a_effect = c_zeros, microbe_b_effect = c_zeros,
  relationship = c_char, relationship_effect = c_zeros,
  therapy = c_char, therapy_effect = c_zeros,
  sun = c_char, medicine = c_char,
  med_green_blue_effect = c_zeros,
  income = c_char, income_effect = c_zeros,
  happiness = c_zeros, total_subjects, study_id_modifier
  ) {

  # Create data frame
  study_id <- paste0(1, study_id_modifier)
  study_id <- c(study_id, study_id, study_id, study_id, study_id)
  df <- data.frame(happiness_original, timepoint,
                   cohort, cohort_effect, study_id,
                   microbe_a, microbe_b, microbe_a_effect,
                   microbe_b_effect, relationship, relationship_effect,
                   therapy, therapy_effect, medicine,
                   med_green_blue_effect, sun,
                   income, income_effect, happiness)

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

  # add income and sun features
  df <- df %>%
    group_by(study_id) %>%
    mutate(income = sample(c("high", "low"), 1),
           sun = sample(c("yes", "no"), 1)) %>%
    ungroup()

  # add relationship and microbe_b features
  df <- df %>%
    group_by(study_id) %>%
    rowwise() %>%
    mutate(relationship = case_when(
           timepoint %in% 3:5 ~ "p",
           timepoint %in% 1:2 ~ sample(c("p", "p", "p", "s"), 1))
      ) %>%
    mutate(microbe_b = case_when(
           timepoint == 1 ~ 0.1,
           timepoint %in% 2:5 ~ 0.5)) %>%
    ungroup()

  # split into thirds; more are affected & in group2
  # for both unaffected and affected
  r <- list_piece(subjects, 3)

  med_unaffected <- r[[1]]
  med_affected <- r[[2]]

  r <- list_piece(med_affected, 3)
  med_affected_1 <- r[[1]]
  med_affected_2 <- r[[2]]

  r <- list_piece(med_unaffected, 3)
  med_unaffected_1 <- r[[1]]
  med_unaffected_2 <- r[[2]]

  # Medication
  df <- df %>%
    group_by(study_id) %>%
    rowwise() %>%
    mutate(medicine = case_when(
      timepoint == 1 ~ sample(c("Pink", "Purple", "Red"), 1),

      # affected individuals (interaction effects)
      study_id %in% med_affected_1 & timepoint %in% 2:3 ~ "Blue",
      study_id %in% med_affected_1 & timepoint %in% 4:5 ~ "Green",

      # THIS IS THE INTERACTION TERM. Only effects happen here (green_blue)
      study_id %in% med_affected_2 & timepoint %in% 2:3 ~ "Green",
      study_id %in% med_affected_2 & timepoint %in% 4:5 ~ "Blue",

      # unaffected individuals (no interaction)
      study_id %in% med_unaffected_1 & timepoint %in% 2:3 ~ "Blue",
      study_id %in% med_unaffected_1 & timepoint %in% 4:5 ~ "Blue",

      study_id %in% med_unaffected_2 & timepoint %in% 2:3 ~ "Green",
      study_id %in% med_unaffected_2 & timepoint %in% 4:5 ~ "Green"
    ))

  add_effects <- function(df, subjects) {

    for (i in subjects) {
      t2_t3 <- subset(df, study_id == i & timepoint %in% 2:3)
      t4_t5 <- subset(df, study_id == i & timepoint %in% 4:5)
      t1 <- subset(df, study_id == i & timepoint == 1)
      t2 <- subset(df, study_id == i & timepoint == 2)
      t3 <- subset(df, study_id == i & timepoint == 3)
      t4 <- subset(df, study_id == i & timepoint == 4)
      t5 <- subset(df, study_id == i & timepoint == 5)

      # medication med_green_blue_effect
      if ("Green" %in% t2_t3$medicine
      && "Blue" %in% t4_t5$medicine) {
        df$med_green_blue_effect[df$timepoint %in% 4:5
        & df$study_id == i] <- my_sample(green_blue_effect, 2)
      }

      # relationship partnered effect
      df$relationship_effect[df$study_id == i
      & df$relationship == "p"] <- my_sample(relationship_p_effect, 1)

      # relationship single effect
      df$relationship_effect[df$study_id == i
      & df$relationship == "s"] <- my_sample(relationship_s_effect, 1)

      # time 1_2 change in relationship
      if ("s" %in% t1$relationship && "p" %in% t2$relationship) {
        df$relationship_effect[df$timepoint == 2
        & df$study_id == i] <- my_sample(relationship_s_p_effect, 2)
      }
      if ("p" %in% t1$relationship && "s" %in% t2$relationship) {
        df$relationship_effect[df$timepoint == 2
        & df$study_id == i] <- my_sample(relationship_p_s_effect, 2)
      }

      # time 1_3 change in relationship
      if ("s" %in% t2$relationship && "p" %in% t3$relationship) {
        df$relationship_effect[df$timepoint == 3
        & df$study_id == i] <- my_sample(relationship_s_p_effect, 2)
      }
      if ("p" %in% t2$relationship && "s" %in% t3$relationship) {
        df$relationship_effect[df$timepoint == 3
        & df$study_id == i] <- my_sample(relationship_p_s_effect, 2)
      }
    }

    # ADD EFFECTS for:
    # cohort, microbe_b, income
    df <- df %>%
      group_by(study_id) %>%
      group_by(timepoint) %>%
      rowwise() %>%
      mutate(cohort_effect = case_when(
             cohort == "A" ~ cohort_a_effect,
             cohort == "B" ~ cohort_b_effect),

             microbe_b_effect = microbe_b,

             income_effect = case_when(
             income == "high" ~ income_high_effect,
             income == "low" ~ income_low_effect)) %>%
      ungroup()

  # effects for microbe_a different in therapy
  # them merge dfs
  df_i1 <- df %>%
    filter(sun == "yes") %>%
    group_by(study_id) %>%
    rowwise() %>%
    mutate(microbe_a_effect = case_when(
           timepoint == 1 ~ my_sample(10, 3),
           timepoint == 2 ~ my_sample(11, 3),
           timepoint == 3 ~ my_sample(20, 3),
           timepoint == 4 ~ my_sample(22, 3),
           timepoint == 5 ~ my_sample(21, 3))) %>%
    mutate(microbe_a_effect = 0.01 * microbe_a_effect,
           microbe_a = microbe_a_effect) %>%
    ungroup()

  df_i2 <- df %>%
    filter(sun == "no") %>%
    group_by(study_id) %>%
    rowwise() %>%
    mutate(microbe_a_effect = 0.1,
           microbe_a = microbe_a_effect) %>%
    ungroup()

  # merge both dataframes back together and sort
  df <- rbind(df_i1, df_i2)
  df <- df %>% arrange(study_id)

  df_i1 <- df %>%
    filter(therapy == "therapy_1") %>%
    group_by(study_id) %>%
    rowwise() %>%
    mutate(therapy_effect = case_when(
           timepoint == 1 ~ my_sample(1, 1),
           timepoint == 2 ~ my_sample(5, 1),
           timepoint == 3 ~ my_sample(10, 1),
           timepoint == 4 ~ my_sample(20, 1),
           timepoint == 5 ~ my_sample(40, 2))) %>%
    ungroup()

  df_i2 <- df %>%
    filter(therapy == "therapy_2") %>%
    group_by(study_id) %>%
    rowwise() %>%
    mutate(therapy_effect = case_when(
           timepoint == 1 ~ my_sample(1, 1),
           timepoint == 2 ~ my_sample(3, 1),
           timepoint == 3 ~ my_sample(5, 1),
           timepoint == 4 ~ my_sample(8, 1),
           timepoint == 5 ~ my_sample(9, 1))) %>%
    ungroup()

  # merge both dataframes back together and sort
  df <- rbind(df_i1, df_i2)
  df <- df %>% arrange(study_id)

  return(df)
  } # end add_effects function

  # remove effect columns and original happiness
  df <- add_effects(df, subjects)
  return(df)
}

# Make cohort A and B dataframes then merge together
df_a <- simulated_df_for_workflow(happiness_original =
rep(happiness_start, total_timepoints),
cohort = rep("A", total_timepoints),
total_subjects = subjects_per_cohort, study_id_modifier = "A")

df_b <- simulated_df_for_workflow(happiness_original =
rep(happiness_start, total_timepoints),
cohort = rep("B", total_timepoints),
total_subjects = subjects_per_cohort, study_id_modifier = "B")

df <- rbind(df_a, df_b)


# Include effect values in happiness/response
df <- df %>%
  mutate(happiness = happiness_original +
  rowSums(.[grep("_effect", names(.))]))

write.table(df, "data/simulation/simulation-effects.txt",
row.names = FALSE, sep = "\t")

# Do not keep effect columns in final dataframe
# the point is to see if vars are selected based on effect
df <- df %>%
  select(-ends_with("_effect")) %>%
  select(!(happiness_original))

# save final dataframe
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

  file_name <- paste0(file_prefix, num_features, "-random-vars-shuffled.txt")
  write.table(df_shuffled_response, file_name, row.names = FALSE, sep = "\t")

  return(df)
}


# add shuffle response option
addFeatures <- function(df, response_var, file_prefix, low, med, hi) {

  df <- read.csv(df, sep = "\t")
  # shuffle response
  shuffled_response <- sample(df[[response_var]])

  # TODO this is repetitive (make list starting with 0)
  df <- writeFiles(file_prefix, df, low, 0, response_var, shuffled_response)
  df <- writeFiles(file_prefix, df, med, low, response_var, shuffled_response)
  df <- writeFiles(file_prefix, df, hi, med, response_var, shuffled_response)

}

addFeatures("data/simulation/simulation.txt", "happiness",
            "data/simulation/simulation-", 10, 100, 1000)
