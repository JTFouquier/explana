
options(warn = -1)
source("scripts/install.R")
package_list <- c("tidyverse", "compositions")
install_r_packages(package_list = package_list)

library(tidyverse)
library(compositions)

ds_name <- snakemake@params[["dataset_name"]]
input_datasets <- snakemake@config[["input_datasets"]]

out_file <- snakemake@output[["out_file"]]
sample_id <- snakemake@config[["sample_id"]]

final_output_folder <- paste0(snakemake@config[["out"]],
                              snakemake@config[["path_dim_preprocess"]],
                              snakemake@params[["dataset_name"]], "/")

method <- snakemake@params[["method"]]
df_mod <- snakemake@params[["df_mod"]]
in_file <- snakemake@input[["in_file"]]

preprocess_dataframe <- function(in_file, method) {

    df <- as.data.frame(readr::read_tsv(in_file))

    if (df_mod != "") {
      eval(parse(text = df_mod))
    }

    if (method %in% c("clr", "arcsin")) {

      s <- df$`sample_id`
      df <- df %>% dplyr::select(-`sample_id`)

      if (method == "clr") {
          df <- as.data.frame(clr(df))
      }

      if (method == "arcsin") {
          arcsin_transform <- function(x) {
              asin(sqrt(x))
          }
          df <- as.data.frame(apply(df, MARGIN = 2, FUN = transform))
      }
      df$`sample_id` <- s
    }
    return(df)
}


check_for_directories <- function(folder_list) {
  for (i in folder_list){
    if (file.exists(i)) {
      break
    } else {
    dir.create(i)
    }
  }
}

preprocess_fp <- paste0(snakemake@config[["out"]],
       snakemake@config[["path_dim_preprocess"]])

check_for_directories(preprocess_fp)

ds_fp <- paste0(snakemake@config[["out"]], final_output_folder)
check_for_directories(ds_fp)

df <- preprocess_dataframe(in_file, method)

write_tsv(df, out_file)
options(warn = 0)
