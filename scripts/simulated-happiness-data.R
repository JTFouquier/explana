# Title     : TODO
# Objective : TODO
# Created by: jenniferfouquier
# Created on: 1/12/22

rm(list = ls())
library("tidyverse")
library("faux")
library("LaplacesDemon")
library(compositions)
require(microbiomeDASim)
library(tidyverse)

out_prefix <- "data/simulation/simulation"

# two cohorts (A and) 
subjects_per_cohort <- 25 # 20 min for testing
happiness_start <- 100 # 100 allows var to show

n_features <- 30
diff_abun_features <- 10

std_dev_random_intercepts <- 7
total_timepoints <- 5 # things might fail if changed

therapy_1_effect <- 40
therapy_2_effect <- 20

# for initializing dataframe
c_zeros <- rep(0, total_timepoints)
c_char <- rep("", total_timepoints)

green_blue_effect <- -20 # -20 works

cohort_a_effect <- -15
cohort_b_effect <- 15

income_high_effect <- 10
income_low_effect <- -10

sun_yes_effect <- 15
sun_no_effect <- -15

s_p_effect <- 12
p_s_effect <- -15

my_sample <- function(x, plus_minus) {
  # get a random number plus or minus input value
  my_min <- x - plus_minus
  my_max <- x + plus_minus
  x <- sample(my_min:my_max, 1)[1]
  return(x)
}


make_microbiome_data <- function(subjects_per_cohort, n_features,
                                 diff_abun_features, sample_ids,
                                 num_timepoints) {
  my_control_mean <- 40
  my_slope <- 30
  df <-
      gen_norm_microbiome(features = n_features, diff_abun_features =
                          diff_abun_features, n_control = subjects_per_cohort,
                          n_treat = subjects_per_cohort, control_mean =
                          my_control_mean, sigma = 15, num_timepoints =  num_timepoints, t_interval = c(1,num_timepoints),
                          rho = 0.7, corr_str = "ar1",func_form = "linear",
                          beta = c(-30, my_slope),missing_pct = 0.0,
                          missing_per_subject = 0, miss_val = NA)

  df <- data.frame(t(df$Y))

  df <- df %>%
    rownames_to_column(var = "Sample")
  write_tsv(df, paste0(out_prefix, "-microbiome-counts.txt"))

  # make compositional (sum to 1)
  df[, -1] <- df[, -1] / rowSums(df[, -1])
  compositional_data <- df[, -1]
  write_tsv(compositional_data,
  paste0(out_prefix, "-microbiome-composition.txt"))

  # clr transform and add IDs (they were ordered)
  clr_transformed <- as.data.frame(clr(compositional_data))
  clr_transformed["sample_id"] <- sample_ids

  write_tsv(clr_transformed, paste0(out_prefix, "-microbiome-clr.txt"))
  return(clr_transformed)
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
  happiness_original = c_zeros, cohort = c_char, cohort_effect = c_zeros,
  timepoint = c(1, 2, 3, 4, 5), study_id = c_char, food_1 = c_zeros,
  food_2 = c_zeros, food_1_effect = c_zeros, food_2_effect = c_zeros,relationship = c_char, relationship_s_p_effect = c_zeros,relationship_p_s_effect = c_zeros, therapy = c_char,
  therapy_1_effect = c_zeros, therapy_2_effect = c_zeros,
  therapy_change_effect = c_zeros, sun = c_char, medicine = c_char,
  med_green_blue_effect = c_zeros,income = c_char, income_effect = c_zeros, happiness = c_zeros, total_subjects, study_id_modifier
  ) {

  # Create data frame
  study_id <- paste0(1, study_id_modifier)
  study_id <- c(study_id, study_id, study_id, study_id, study_id)
  df <- data.frame(happiness_original, timepoint, cohort, cohort_effect,study_id, food_1, food_2, food_1_effect, food_2_effect, relationship,
                   relationship_s_p_effect, relationship_p_s_effect,
                   therapy, therapy_change_effect, therapy_1_effect,
                   therapy_2_effect, medicine, med_green_blue_effect, sun,
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

  df_cohort_a <- base_dataframe(df, cohort_a_subjects, subjects, "cohort_a")
  df_cohort_b <- base_dataframe(df, cohort_b_subjects, subjects, "cohort_b")
  df <- rbind(df_cohort_a, df_cohort_b)

  # add income and sun features
  df <- df %>%
    group_by(study_id) %>%
    mutate(income = sample(c("high", "low"), 1),
           sun = sample(c("yes", "yes", "yes", "yes", "no", "no"), 1)) %>%
    ungroup()

  # add relationship and food_2 features
  df <- df %>%
    group_by(study_id) %>%
    rowwise() %>%
    mutate(relationship = case_when(
           timepoint %in% 3:5 ~ "p",
           timepoint %in% 1:2 ~ sample(c("p", "p", "p", "s", "s"), 1))
      ) %>%
    mutate(food_2 = case_when(
           timepoint == 1 ~ 1,
           timepoint %in% 2:5 ~ 5)) %>%
    ungroup()

  # affected means switched pill color; otherwise the same
  # split by factor x list_piece function
  r <- list_piece(subjects, 3)

  med_unaffected <- r[[1]]
  med_affected <- r[[2]]

  r <- list_piece(med_affected, 3)
  med_affected_1 <- r[[1]]
  med_affected_2 <- r[[2]]  # green_blue effect group

  r <- list_piece(med_unaffected, 2)
  med_unaffected_1 <- r[[1]]
  med_unaffected_2 <- r[[2]]

  # Medication
  df <- df %>%
    group_by(study_id) %>%
    rowwise() %>%
    mutate(medicine = case_when(
      timepoint == 1 ~ "blue",

      study_id %in% med_affected_1 & timepoint %in% 2:3 ~ "blue",
      study_id %in% med_affected_1 & timepoint %in% 4:5 ~ "green",

      # THIS IS THE INTERACTION TERM. Only effects happen here (green_blue)
      study_id %in% med_affected_2 & timepoint %in% 2:3 ~ "green",
      study_id %in% med_affected_2 & timepoint %in% 4:5 ~ "blue",

      study_id %in% med_unaffected_1 & timepoint %in% 2:3 ~ "blue",
      study_id %in% med_unaffected_1 & timepoint %in% 4:5 ~ "blue",

      study_id %in% med_unaffected_2 & timepoint %in% 2:3 ~ "green",
      study_id %in% med_unaffected_2 & timepoint %in% 4:5 ~ "green"
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
      if ("green" %in% t2_t3$medicine
      && "blue" %in% t4_t5$medicine) {
        df$med_green_blue_effect[df$timepoint %in% 4:5
        & df$study_id == i] <- my_sample(green_blue_effect, 2)
      }

      # time 1_2 change in relationship
      if ("s" %in% t1$relationship && "p" %in% t2$relationship) {
        df$relationship_s_p_effect[df$timepoint == 2
        & df$study_id == i] <- my_sample(s_p_effect, 1)
      }
      if ("p" %in% t1$relationship && "s" %in% t2$relationship) {
        df$relationship_p_s_effect[df$timepoint == 2
        & df$study_id == i] <- my_sample(p_s_effect, 1)
      }

      # time 1_3 change in relationship
      if ("s" %in% t2$relationship && "p" %in% t3$relationship) {
        df$relationship_s_p_effect[df$timepoint == 3
        & df$study_id == i] <- my_sample(s_p_effect, 1)
      }
      if ("p" %in% t2$relationship && "s" %in% t3$relationship) {
        df$relationship_p_s_effect[df$timepoint == 3
        & df$study_id == i] <- my_sample(p_s_effect, 1)
      }
    }

    # ADD EFFECTS for:
    # cohort, food_2, income
    df <- df %>%
      group_by(study_id) %>%
      group_by(timepoint) %>%
      rowwise() %>%
      mutate(cohort_effect = case_when(
             cohort == "cohort_a" ~ cohort_a_effect,
             cohort == "cohort_b" ~ cohort_b_effect),

             sun_effect = case_when(
             sun == "yes" ~ sun_yes_effect,
             sun == "no" ~ sun_no_effect),

             food_1 = case_when(
             timepoint == 1 ~ my_sample(10, 2),
             timepoint == 2 ~ my_sample(25, 2),
             timepoint == 3 ~ my_sample(40, 2),
             timepoint == 4 ~ my_sample(55, 2),
             timepoint == 5 ~ my_sample(70, 2)),

             food_1_effect = food_1 * 0.01,

             # TODO try to do this in only some people
             food_2_effect = case_when(
             timepoint == 2 ~ 12,
             timepoint != 2 ~ 1),

             income_effect = case_when(
             income == "high" ~ income_high_effect,
             income == "low" ~ income_low_effect)) %>%
      ungroup()

  # effects are stronger in therapy_1
  df_i1 <- df %>%
    filter(therapy == "therapy_1") %>%
    group_by(study_id) %>%
    rowwise() %>%
    mutate(therapy_1_effect = case_when(
           timepoint == 1 ~ 0,
           timepoint != 1 ~ therapy_1_effect),
           therapy_change_effect = case_when(
           timepoint == 1 ~ my_sample(30, 3),
           timepoint == 2 ~ my_sample(60, 3),
           timepoint == 3 ~ my_sample(90, 3),
           timepoint == 4 ~ my_sample(120, 3),
           timepoint == 5 ~ my_sample(150, 3))) %>%
    ungroup()

  df_i2 <- df %>%
    filter(therapy == "therapy_2") %>%
    group_by(study_id) %>%
    rowwise() %>%
    mutate(therapy_2_effect = case_when(
           timepoint == 1 ~ 0,
           timepoint != 1 ~ therapy_2_effect),
           therapy_change_effect = case_when(
           timepoint == 1 ~ my_sample(10, 4),
           timepoint == 2 ~ my_sample(15, 4),
           timepoint == 3 ~ my_sample(20, 4),
           timepoint == 4 ~ my_sample(25, 4),
           timepoint == 5 ~ my_sample(30, 4))) %>%
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


make_random_df <- function(df, distribution_names, num_cols,
                           cols_in_random_df) {

  n_rows <- subjects_per_cohort * 2 * total_timepoints
  # code for data distribution obtained using ChatGPT
  set.seed(111)
  data_list <- list()
  for (distribution in distribution_names) {
    for (i in 1:num_cols) {
      column_name <- paste0("V_", distribution, "_", i)
      if (distribution == "normal") {
        data_list[[column_name]] <- rnorm(n_rows)
      } else if (distribution == "bernoulli") {
        data_list[[column_name]] <- rbinom(n_rows, 1, 0.5)
      } else if (distribution == "binomial") {
        data_list[[column_name]] <- rbinom(n_rows, 10, 0.3)
      } else if (distribution == "poisson") {
        data_list[[column_name]] <- rpois(n_rows, 2)
      } else if (distribution == "exponential") {
        data_list[[column_name]] <- rexp(n_rows, 0.5)
      } else if (distribution == "gamma") {
        data_list[[column_name]] <- rgamma(n_rows, 2, 1)
      } else if (distribution == "weibull") {
        data_list[[column_name]] <- rweibull(n_rows, 2, 1)
      } else if (distribution == "dirichlet") {
        data_list[[column_name]] <- rdirichlet(n_rows,
        rep(1, cols_in_random_df))  # total "microbes/compositional"
      } else {
        stop(paste0("distribution unavailable: ", distribution))
      }
    }
  }
  # generate the random df, give name and merge
  random_df <- as.data.frame(data_list)
  random_df[["sample_id"]] <- df[["sample_id"]]
  df <- merge(df, random_df, by.y = "sample_id", by.x = "sample_id")
  return(df)
}


writeFiles <- function(file_prefix, df, n_features, n_dirichlet,
response_var, sampled) {

  # Dirichlet distribution on its own
  dirichlet_distribution <- c("dirichlet")
  cols_in_random_df <- n_dirichlet  # this is how many "microbes"
  num_cols <- ceiling((cols_in_random_df / cols_in_random_df))
  df <- make_random_df(df, dirichlet_distribution, num_cols, cols_in_random_df)

  # Other data distributions
  other_distributions <- c("normal", "bernoulli", "binomial",
  "poisson", "exponential", "gamma", "weibull")

  cols_in_random_df <- length(other_distributions)
  num_cols <- ceiling((n_features / cols_in_random_df))
  df <- make_random_df(df, other_distributions, num_cols, cols_in_random_df)

  df_shuffled_response <- df
  df_shuffled_response[[response_var]] <- sampled

  # get actual total number of features
  n_features_all_distributions <-
  (num_cols * length(other_distributions)) + n_dirichlet

  file_name <- paste0(file_prefix,
  n_features_all_distributions, "-random-vars.txt")
  write.table(df, file_name, row.names = FALSE, sep = "\t")

  file_name <- paste0(file_prefix,
  n_features_all_distributions, "-random-vars-shuffled.txt")
  write.table(df_shuffled_response, file_name, row.names = FALSE, sep = "\t")
}


# add shuffle response option
addFeatures <- function(df, response_var, file_prefix,
feature_count_list, n_dirichlet) {

  # df <- read.csv(df, sep = "\t")
  # shuffle response
  shuffled_response <- sample(df[[response_var]])

  for (i in feature_count_list) {
    writeFiles(file_prefix, df, i, n_dirichlet,
    response_var, shuffled_response)
  }
}

# Make cohort A and B dataframes then merge together

df_a <- simulated_df_for_workflow(happiness_original =
rep(happiness_start, total_timepoints),
cohort = rep("cohort_a", total_timepoints),
total_subjects = subjects_per_cohort, study_id_modifier = "a")

df_b <- simulated_df_for_workflow(happiness_original =
rep(happiness_start, total_timepoints),
cohort = rep("cohort_b", total_timepoints),
total_subjects = subjects_per_cohort, study_id_modifier = "b")

df <- rbind(df_a, df_b)

# Include effect values in happiness/response
df <- df %>%
  mutate(happiness = happiness_original +
  rowSums(.[grep("_effect", names(.))]))

write.table(df, paste0(out_prefix, "-effects.txt"),
row.names = FALSE, sep = "\t")

# Do not keep effect columns in final dataframe
# the point is to see if vars are selected based on effect
df <- df %>%
  select(-ends_with("_effect")) %>%
  select(!(happiness_original))

# sort by therapy for merging with simulated data; save df
df <- df %>%
  arrange(desc(therapy), sample_id, timepoint)
write.table(df, paste0(out_prefix, ".txt"), row.names = FALSE, sep = "\t")

df_microbes <- make_microbiome_data(subjects_per_cohort, n_features,
diff_abun_features, df["sample_id"], total_timepoints)

df_microbes <- merge(df, df_microbes)

df_microbes <- df_microbes %>%
 arrange(desc(therapy), sample_id, timepoint)
write.table(df_microbes, paste0(out_prefix, "-microbiome.txt"),
row.names = FALSE, sep = "\t")

addFeatures(df_microbes, "happiness", out_prefix,
            feature_count_list = c(105, 315),  # feature n
            n_dirichlet = 50)



