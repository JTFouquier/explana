library(tidyverse)
library(MASS)
library(ape)
library(purrr)
library(ggplot2)
library(rstatix)
library(vars)

color_list <- list('First' = 'cornflowerblue',
                    'Previous' = 'lightgreen',
                    'Pairwise' = '#AB82FF',
                    'Original' = '#FFFACD')

postHocGraphs = function(metadata_long, gg_title){

  print(metadata_long)
  g = ggplot(metadata_long, aes(x = PredictorValue, y = get(response))) +
  ylab(response) +
  theme(strip.background = element_rect(fill=color_list[[gg_info]]))
  return(g)
}

splitMetadata = function(plotting_vars, df_wide, cols_for_graph){
  # important features are the ones to plot
  informative_vars = c(cols_for_graph, plotting_vars)

  # only keep necessary columns; long format needed for facet wrap
  df_small = df_wide[informative_vars]

  df_long = gather(df_small, key = "Predictor",
                   value = "PredictorValue", all_of(plotting_vars))
  return(df_long)
}


df = read.csv(df, sep = "\t", check.names=FALSE)
# TODO fix timepoint
df$Timepoint = as.factor(df$Timepoint)

important_features = read.csv(important_features, sep = "\t")
important_features = unique(important_features$decoded.features)

# These are not part of the facets -- just needed for display, etc
cols_for_graph = c(random_effect, response, "Timepoint")
important_columns = c(cols_for_graph, important_features)
facet.cols = 2 # TODO fix hardcoding here

# get dataframes with important vars and either numeric or categoric predictors and convert to long
df_col_selection = df[important_features]


df_categorical = df_col_selection %>% select(where(negate(is.numeric)))
df_numerical = df_col_selection %>% select(where(is.numeric))

categorical_vars = colnames(df_categorical)
numerical_vars = colnames(df_numerical)

if (length(categorical_vars) > 0){
  df_long_categorical = splitMetadata(categorical_vars, df, cols_for_graph)
  g = postHocGraphs(df_long_categorical, gg_title)
  g = g +
   geom_point(position=position_jitter(width=0.1, height=0.1)) +
   geom_boxplot(outlier.colour=NA, fill=NA, colour="grey60") +
   facet_wrap(vars(Predictor), ncol = facet.cols, scales = "free_x") +
   labs(title = paste0(gg_info, " (Categorical Variables)"))
print(g)
ggsave(out_file_categorical, g)
} else {
  return()
}

if (length(numerical_vars) > 0){
  df_long_numerical = splitMetadata(numerical_vars, df, cols_for_graph)
  g = postHocGraphs(df_long_numerical, gg_title)
  g = g +
    geom_point() +
    geom_smooth(method='lm', se = FALSE) +
    facet_wrap(vars(Predictor), ncol = facet.cols, scales = "free_x") +
    labs(title = paste0(gg_info, " (Numerical Variables)"))
  print(g)
  ggsave(out_file, g)
} else {
  return()
}



# expand color palettes
#color.count = length(unique(df$StudyID))
#mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(color.count)




