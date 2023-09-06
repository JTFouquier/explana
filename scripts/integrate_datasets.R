source("scripts/install.R")
package_list <- c("dplyr")
install_r_packages(package_list = package_list)

library(dplyr)

source("scripts/viz-datatable.R")


output_folder <- paste0(snakemake@config[["out"]],
snakemake@config[["path_original"]])

sample_id <- snakemake@config[["sample_id"]]
study_id <- snakemake@config[["random_effect"]]
response_var <- snakemake@config[["response_var"]]

file_path_list <- snakemake@input[["file_path_list"]]

build_datatable <- snakemake@params[["build_datatable"]]
ds_param_dict_list <- snakemake@params[["ds_param_dict_list"]]
timepoint_config <- snakemake@config[["timepoint"]]


verify_inputs <- function(constrain_cols, drop_cols) {
    if (length(constrain_cols) > 0 & length(drop_cols) > 0) {
    # If both these fields exist the entire script fails; TODO make error
    print("WORKFLOW WARNING: Cannot have features (arguments) in 'drop_cols'
    and 'constrain_cols' parameters")
    return(FALSE)
  }
  return(TRUE)
}


check_duplicate_colnames <- function(df, df_file_name) {
  complete_column_list <- list(colnames(df))
  complete_column_list[sample_id] <- NULL
  if (length(complete_column_list) != length(unique(complete_column_list))) {
    print(paste0("WORKFLOW WARNING: duplicate column names found in dataset ",
                 df_file_name, " will affect analysis."))
  }
}


filter_dataframe <- function(df, df_mod_list, df_complete_flag) {

  # For complete/full/merged dataframe, get modification information from
  # Config file. For individual datasets, get from json
  if (df_complete_flag == TRUE) {
    # for after datasets are merged together
    drop_rows <- eval(parse(text = snakemake@config[["drop_rows"]]))
    constrain_rows <- eval(parse(text = snakemake@config[["constrain_rows"]]))
    drop_cols <- eval(parse(text = snakemake@config[["drop_cols"]]))
    constrain_cols <- eval(parse(text = snakemake@config[["constrain_cols"]]))
  } else if (df_complete_flag == FALSE) {
    # for individual datasets
    constrain_cols <- df_mod_list[["constrain_cols"]]
    drop_cols <- df_mod_list[["drop_cols"]]
    constrain_rows <- df_mod_list[["constrain_rows"]]
    drop_rows <- df_mod_list[["drop_rows"]]
  }
  verify_flag <- verify_inputs(constrain_cols = constrain_cols,
                              drop_cols = drop_cols)
  if (verify_flag == FALSE) {
    return()
  }

  check_duplicate_colnames(df, df_file_name)

  if (length(constrain_cols) > 0) {
    df <- df %>% select(all_of(constrain_cols))
  }
  if (length(drop_cols) > 0) {
    df <- df %>% select(!all_of(drop_cols))
  }
  if (length(drop_rows) > 0) {
    for (i in names(drop_rows)) {
      df <- df %>% filter(.data[[{{i}}]] != drop_rows[[i]])
    }
  }
  if (length(constrain_rows) > 0) {
    for (i in names(constrain_rows)) {
      df <- df %>% filter(.data[[{{i}}]] == constrain_rows[[i]])
    }
  }
  return(df)
}


main <- function(file_path_list, ds_param_dict_list) {
  if (length(file_path_list) != length(ds_param_dict_list)) {
    print(print(paste0("WORKFLOW WARNING: Dataset count does not equal
    dataset parameter dictionary count")))
  }
  df_list <- list()
  for (i in 1:length(file_path_list)) {
    # For each dataset, get the parameter dict (drop, constrain, etc)
    ds_param_dict <- ds_param_dict_list[[i]]
    df_mod_list <- snakemake@config[[ds_param_dict]]
    df_mod_list <- eval(parse(text = df_mod_list))
    df_complete_flag <- FALSE
    df <- as.data.frame(read_tsv(file_path_list[[i]]))
    df <- filter_dataframe(df, df_mod_list, df_complete_flag)
    df_list <- append(df_list, list(df))
  }
  df_complete_flag <- TRUE
  df_complete <- df_list %>% reduce(full_join, by = sample_id)
  df_complete <- filter_dataframe(df_complete, df_mod_list, df_complete_flag)

  # make sure the timepoint column is ranked (ordered)
  df_complete$timepoint_original <- df_complete[[timepoint_config]]

  df_complete <- df_complete %>%
    arrange(timepoint_original) %>%
    mutate(timepoint_explana = dense_rank(timepoint_original)) %>%
    select(-{{ timepoint_config }})

  write_tsv(df_complete, paste0(output_folder, "original.txt"))
}


main(file_path_list, ds_param_dict_list)

# TODO optional create visualizer
# TODO change this to either df or load the dataframe
if (build_datatable %in% c("TRUE", "True")) {
  input_file_name <- paste0(output_folder, "original.txt")
  output_file_name <- paste0(output_folder, "vizualizer-merged-dfs.html")
  build_datatable_viz(input_file_name, output_file_name)
}
