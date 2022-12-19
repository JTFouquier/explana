

library(tidyverse)
source("scripts/viz-datatable.R")

output_folder = paste0(snakemake@config[["out"]], snakemake@config[["path_merged_data"]])
sample_id = snakemake@config[["sample_id"]]
file_path_list = snakemake@input[["file_path_list"]]
build_datatable = snakemake@params[["build_datatable"]]

drop_rows = eval(parse(text=snakemake@config[["drop_rows"]]))
constrain_rows = eval(parse(text=snakemake@config[["constrain_rows"]]))
drop_cols = eval(parse(text=snakemake@config[["drop_cols"]]))
constrain_cols = eval(parse(text=snakemake@config[["constrain_cols"]]))


filter_dataframe = function(){
    if (length(constrain_cols) > 0 & length(drop_cols) > 0){
    print("Cannot have features (arguments) in 'drop_cols' and 'constrain_cols'
    parameters")
    return()
  }

  unique_cols = list()
  for (df_fp in file_path_list){
      df = read.table(df_fp,sep="\t",header=T)
      unique_cols = append(unique_cols, colnames(df))
  }
  unique_cols["SampleID"] = NULL
  if (length(unique_cols) != length(unique(unique_cols))){
    print("Warning, duplicate column names in datasets could affect analysis")
  }

  df_list <- lapply(file_path_list, read_tsv)
  df_large = df_list %>% reduce(full_join, by='SampleID')

  if (length(constrain_cols) > 0){
    df_large = df_large %>%
      select(all_of(constrain_cols))
  }
  if (length(drop_cols) > 0){
    df_large = df_large %>%
      select(!all_of(drop_cols))
  }
  for (i in names(drop_rows)){
    df_large = df_large %>%
        filter(.data[[{{i}}]] != drop_rows[[i]])
  }
  for (i in names(constrain_rows)){
    df_large = df_large %>%
        filter(.data[[{{i}}]] == constrain_rows[[i]])
  }
  write_tsv(df_large, paste0(output_folder, "final-merged-dfs.txt"))
  return()
}


filter_dataframe()

# TODO optional create visualizer
# TODO change this to either df or load the dataframe
if (build_datatable %in% c("TRUE", "True")){
  input_file_name = paste0(output_folder, "final-merged-dfs.txt")
  output_file_name = paste0(output_folder, "vizualizer-merged-dfs.html")
  build_datatable_viz(input_file_name, output_file_name)
}
