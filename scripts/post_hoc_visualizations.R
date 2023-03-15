package_list <- c("tidyverse", "MASS", "ape", "purrr", "ggplot2",
"rstatix", "vars")

package_list <- package_list[!(package_list %in%
installed.packages()[, "Package"])]
if (length(package_list) > 0) {
  install.packages(package_list)
}

library(tidyverse)
library(MASS)
library(ape)
library(purrr)
library(ggplot2)
library(rstatix)
library(vars)

color_list <- list("First" = "cornflowerblue",
                   "Previous" = "lightgreen",
                   "Pairwise" = "#AB82FF",
                   "Original" = "#FFFACD")


wide_to_long <- function(plotting_vars, df_wide, cols_for_graph) {
  # important features are the ones to plot
  informative_vars <- c(cols_for_graph, plotting_vars)

  # only keep necessary columns; long format needed for facet wrap
  df_small <- df_wide[informative_vars]

  df_long <- gather(df_small, key = "Predictor",
                    value = "PredictorValue", all_of(plotting_vars))
  return(df_long)
}


# set up graphs for both numerical and categorical graphs
post_hoc_graphs <- function(metadata_long, gg_info) {
  g <- ggplot(metadata_long, aes(x = PredictorValue, y = get(response))) +
  ylab(response) +
  theme(strip.background = element_rect(fill = color_list[[gg_info]]))
  return(g)
}


plot_cat_vars <- function(categorical_vars, df, cols_for_graph, gg_info,
                          df_long_categorical, facet_cols) {
  if (length(categorical_vars) > 0) {
    df_long_categorical <- wide_to_long(categorical_vars, df, cols_for_graph)
    g <- post_hoc_graphs(df_long_categorical, gg_info)
    g <- g +
    geom_point(position = position_jitter(width = 0.1, height = 0.1)) +
    geom_boxplot(outlier.colour = NA, fill = NA, colour = "grey60") +
    facet_wrap(vars(Predictor), ncol = facet_cols, scales = "free_x") +
    labs(title = paste0(gg_info, " (Categorical Variables)"))
  return(g)
  } else {
    return()
  }
}


plot_num_vars <- function(numerical_vars, df, cols_for_graph, gg_info,
                          df_long_numerical, facet_cols) {
  if (length(numerical_vars) > 0) {
    df_long_numerical <- wide_to_long(numerical_vars, df, cols_for_graph)
    g <- post_hoc_graphs(df_long_numerical, gg_info)
    g <- g +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      facet_wrap(vars(Predictor), ncol = facet_cols, scales = "free_x") +
      labs(title = paste0(gg_info, " (Numerical Variables)"))
    return(g)
  } else {
    return()
  }
}


post_hoc_viz <- function(df, important_features, gg_info, df_path) {

  out_file <-
    paste0(snakemake@config[["out"]], df_path, "post-hoc-visualizations.pdf")

  df <- read.csv(df, sep = "\t", check.names = FALSE)
  # TODO fix timepoint
  df$Timepoint <- as.factor(df$Timepoint)

  important_features <- read.csv(important_features, sep = "\t")
  important_features <- unique(important_features$decoded.features)

  # These are not part of the facets -- just needed for display, etc
  cols_for_graph <- c(random_effect, response, "Timepoint")
  important_columns <- c(cols_for_graph, important_features)
  facet_cols <- 2 # TODO fix hardcoding here

  # get dataframes with important vars and either numeric or categoric predictors and convert to long
  df_col_selection <- df[important_features]

  df_categorical <- df_col_selection %>% select(where(negate(is.numeric)))
  df_numerical <- df_col_selection %>% select(where(is.numeric))

  categorical_vars <- colnames(df_categorical)
  numerical_vars <- colnames(df_numerical)

  g_cat <- plot_cat_vars(categorical_vars, df, cols_for_graph, gg_info,
  df_long_categorical, facet_cols)
  g_num <- plot_num_vars(numerical_vars, df, cols_for_graph, gg_info,
  df_long_numerical, facet_cols)

  g_complete <- ggarrange(g_cat, g_num, ncol = 2, nrow = 1)

  ggsave(out_file, g_complete, width = 8.5, height = 11)
  return(out_file)
}



# expand color palettes
#color.count = length(unique(df$StudyID))
#mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(color.count)