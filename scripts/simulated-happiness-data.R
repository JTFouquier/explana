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
out_prefix <- "data/simulation/simulation-revert"

# two cohorts (A and)
subjects_per_cohort <- 30 # 20 min for testing
happiness_start <- 100 # 100 allows var to show

# simulated longitudinal microbiome using microbiomeDASim
n_features <- 30
diff_abun_features <- 5

total_timepoints <- 5 # things might fail if changed

# for initializing dataframe
c_zeros <- rep(0, total_timepoints)
c_char <- rep("", total_timepoints)

green_blue_effect <- -22 # -20 works

unhealthy_effect <- -10
healthy_effect <- 11

income_high_effect <- 10
income_low_effect <- -11

vacation_yes_effect <- 10
vacation_no_effect <- -11

s_p_effect <- 7
p_s_effect <- -13

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
                          n_treat = subjects_per_cohort,
                          control_mean = my_control_mean, sigma = 5,
                          num_timepoints = num_timepoints,
                          t_interval = c(1, num_timepoints),
                          rho = 0.7, corr_str = "ar1", func_form = "linear",
                          beta = c(-30, my_slope), missing_pct = 0.0,
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
    mutate(intervention = case_when(
      study_id %in% subjects_half1 ~ "treated",
      study_id %in% subjects_half2 ~ "control"
    ))

  return(df)
}


add_random_effects_subjects <- function(df) {

  between <- list(intervention = c(treated = "treated",
                                   control = "control"))

  within <- list(timepoint = c(1, 2, 3, 4, 5))
  mu <- data.frame(treated = c(35, 65, 95, 125, 155),
                   control = c(35, 35, 35, 35, 35),
                   row.names = within$timepoint)
  # simulated dataset
  df_sim <- sim_design(within, between,
                       n = subjects_per_cohort, mu = mu,
                       sd = 5, r = 0.7, plot = FALSE)

  df_sim <- df_sim %>%
    pivot_longer(cols = -c(id, intervention),
    names_to = "timepoint", values_to = "intervention_change_effect") %>%
    mutate(timepoint = as.integer(timepoint))

  df_sim <- df_sim %>%
    arrange(intervention) %>%
    select(intervention_change_effect)

  df <- df %>%
    arrange(desc(intervention)) %>%
    select(!intervention_change_effect)

  df["intervention_change_effect"] <- df_sim["intervention_change_effect"]
  return(df)
}


# Some are placeholders here, but modified later
simulated_df_for_workflow <- function(cohort = c_char,
                                      total_subjects,
                                      study_id_modifier) {

  happiness_original <- rep(happiness_start, total_timepoints)

  # Create data frame
  study_id <- paste0(1, study_id_modifier)
  study_id <- c(study_id, study_id, study_id, study_id, study_id)


  df <- data.frame(happiness_original, cohort, cohort_effect = c_zeros,
                   timepoint = c(1, 2, 3, 4, 5), study_id, sunshine = c_zeros,
                   relaxation = c_zeros, sunshine_effect = c_zeros,
                   relaxation_effect = c_zeros, relationship = c_char,
                   relationship_s_p_effect = c_zeros,
                   relationship_p_s_effect = c_zeros, intervention = c_char,
                   intervention_change_effect = c_zeros, vacation = c_char,
                   vacation_effect = c_zeros,
                   medicine = c_char, med_green_blue_effect = c_zeros,
                   income = c_char, income_effect = c_zeros,
                   happiness = c_zeros)

  num_columns <- length(colnames(df))

  new_subject_df <- df
  for (i in 2:total_subjects){
      new_subject_df$study_id <- paste0(i, study_id_modifier)
      df <- rbind(df, new_subject_df[rep(1:num_columns, 1)])
  }

  subjects <- unique(df$study_id)

  df$sample_id <- paste0(df$study_id, ".", df$timepoint)
  unhealthy_subjects <- subjects[1:(length(subjects) / 2)]
  healthy_subjects <- subjects[!subjects %in% unhealthy_subjects]

  df_unhealthy <- base_dataframe(df, unhealthy_subjects, subjects, "unhealthy")
  df_healthy <- base_dataframe(df, healthy_subjects, subjects, "healthy")
  df <- rbind(df_unhealthy, df_healthy)

  df <- df %>%
    group_by(study_id) %>%
    mutate(income = sample(c("high", "low"), 1),
           vacation = sample(c("yes", "yes", "no"), 1)) %>%
    ungroup()

  df <- df %>%
    group_by(study_id) %>%
    rowwise() %>%
    mutate(relationship = case_when(
           timepoint %in% 3:5 ~ "p",
           timepoint %in% 1:2 ~ sample(c("p", "p", "p", "s", "s"), 1))
      ) %>%
    mutate(relaxation = case_when(
           timepoint == 1 ~ 2,
           timepoint %in% 2:5 ~ 12)) %>%
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
    # cohort, relaxation, income
    df <- df %>%
      group_by(study_id) %>%
      group_by(timepoint) %>%
      rowwise() %>%
      mutate(cohort_effect = case_when(
             cohort == "unhealthy" ~ unhealthy_effect,
             cohort == "healthy" ~ healthy_effect),

             vacation_effect = case_when(
             vacation == "yes" ~ vacation_yes_effect,
             vacation == "no" ~ vacation_no_effect),

             sunshine = case_when(
             timepoint == 1 ~ 10,
             timepoint == 2 ~ 25,
             timepoint == 3 ~ 40,
             timepoint == 4 ~ 55,
             timepoint == 5 ~ 70),

             sunshine_effect = sunshine * 0.1,
             relaxation_effect = relaxation,

             income_effect = case_when(
             income == "high" ~ income_high_effect,
             income == "low" ~ income_low_effect)) %>%
      ungroup()

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
add_noisy_variables <- function(df, response_var, file_prefix,
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

df_a <- simulated_df_for_workflow(cohort = rep("unhealthy", total_timepoints),
total_subjects = subjects_per_cohort, study_id_modifier = "unhealthy")

df_b <- simulated_df_for_workflow(cohort = rep("healthy", total_timepoints),
total_subjects = subjects_per_cohort, study_id_modifier = "healthy")

df <- rbind(df_a, df_b)

# use package faux to add mixed effects to simulation
df <- add_random_effects_subjects(df)

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

# sort by intervention for merging with simulated data; save df
df <- df %>%
  arrange(intervention, sample_id, timepoint)
write.table(df, paste0(out_prefix, ".txt"), row.names = FALSE, sep = "\t")

df_microbes <- make_microbiome_data(subjects_per_cohort, n_features,
diff_abun_features, df["sample_id"], total_timepoints)

df_microbes <- merge(df, df_microbes)

df_microbes <- df_microbes %>%
 arrange(intervention, sample_id, timepoint)
write.table(df_microbes, paste0(out_prefix, "-microbiome.txt"),
row.names = FALSE, sep = "\t")

add_noisy_variables(df_microbes, "happiness", out_prefix,
            feature_count_list = c(105, 405),  # feature n
            n_dirichlet = 95)
