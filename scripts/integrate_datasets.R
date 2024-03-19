source("scripts/install.R")
package_list <- c("dplyr", "readr")
install_r_packages(package_list = package_list)

library(readr)
library(dplyr)


output_folder <- paste0(snakemake@config[["out"]],
snakemake@config[["path_original"]])

sample_id <- snakemake@config[["sample_id"]]
study_id <- snakemake@config[["random_effect"]]
response_var <- snakemake@config[["response_var"]]

file_path_list <- snakemake@input[["file_path_list"]]

build_datatable <- snakemake@params[["build_datatable"]]
ds_param_dict_list <- snakemake@params[["ds_param_dict_list"]]
timepoint_config <- snakemake@config[["timepoint"]]


check_duplicate_colnames <- function(df, df_file_name) {
  complete_column_list <- list(colnames(df))
  complete_column_list[sample_id] <- NULL
  if (length(complete_column_list) != length(unique(complete_column_list))) {
    print(paste0("WORKFLOW WARNING: duplicate column names found in dataset ",
                 df_file_name, " will affect analysis."))
  }
}


filter_dataframe <- function(df, df_mod_list) {
  df_mod <- snakemake@config[["df_mod"]]
  eval(parse(text = df_mod))
  check_duplicate_colnames(df, df_file_name)
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
    df_mod_list <- ds_param_dict
    df <- as.data.frame(readr::read_tsv(file_path_list[[i]]))
    df_list <- append(df_list, list(df))
  }

  df_complete <- df_list %>% purrr::reduce(inner_join, by = sample_id)
  df_complete <- filter_dataframe(df_complete, df_mod_list)

  # make sure the timepoint column is ranked (ordered)
  df_complete$timepoint_original <- df_complete[[timepoint_config]]

  df_complete <- df_complete %>%
    arrange(timepoint_original) %>%
    mutate(timepoint_explana = dense_rank(timepoint_original)) %>%
    select(-{{ timepoint_config }})

  write_tsv(df_complete, paste0(output_folder, "original.txt"))
}


main(file_path_list, ds_param_dict_list)
