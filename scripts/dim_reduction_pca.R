# Script that assists with performing dimensionality reduction via
# principal component analysis (PCA) on sets of user specified input variables

# Important PCs are selected past cumulative threshold cutoff.
# The most important PCs replace the original variables in the dataframe.
#
# var1 var2 var3 might get converted to dim1 and dim2;
# Overall, a dimensionality reduction from 3 to 2 variables
source("scripts/install.R")
package_list <- c("factoextra", "MASS", "dplyr")
install_r_packages(package_list = package_list)

library(factoextra)
library(MASS)
library(dplyr)

ds_name <- snakemake@params[["dataset_name"]]
input_datasets <- snakemake@config[["input_datasets"]]
in_file <- input_datasets[[ds_name]][["file_path"]]

out_file <- snakemake@output[["out_file"]]
sample_id <- snakemake@config[["sample_id"]]

# Need to get the name of the user declared 'pca_groups_list' using method
# from Snakefile, then get the actual data object from config file.
pca_list_name <- snakemake@params[["pca_list"]]
pca_list <- pca_list_name

final_output_folder <- paste0(snakemake@config[["out"]],
                              snakemake@config[["path_dim_pca"]],
                              snakemake@params[["dataset_name"]], "/")


write_table_special <- function(df, folder_name, file_name) {
  outf <- paste0(folder_name, file_name)
  write.table(df, outf, row.names = FALSE, sep = "\t")
}
inf <- paste0(in_file)
df <- read.table(inf, header = TRUE, sep = "\t", strip.white = FALSE)

check_for_directories <- function(folder_list) {
  for (i in folder_list){
    if (file.exists(i)) {
      break
    } else {
    dir.create(i)
    }
  }
}


perform_pca <- function(pca_list, pca_name, df, final_output_folder) {

  # get items from each list describing how to perform pca
  pca_list <- list(pca_list[pca_name][[1]])[[1]]

  output_folder <- paste0(final_output_folder, pca_name, "/")
  feature_list <- eval(parse(text = pca_list$feature_list))
  cum_var_thres <- eval(parse(text = pca_list$cum_var_thres))

  pca_df <- df %>% select(all_of(feature_list))
  # Scale PCA plot
  # TODO prcomp. this brings up need to allow user adjustments to functions
  pca_scaled <- prcomp(pca_df, scale = TRUE)

  eig_val <- get_eigenvalue(pca_scaled)
  eig_val # data frame with all info on eigen values and variance of axes

  res_var <- get_pca_var(pca_scaled)
  res_var # data frame of data frames with all the info on the variables
  var_coord <- as.data.frame(res_var$coord)     # Coordinates
  var_cor <- as.data.frame(res_var$cor)         # Corrs between vars and dims
  var_contrib <- as.data.frame(res_var$contrib) # Contributions to the PCs
  var_cos2 <- as.data.frame(res_var$cos2)       # Quality of representation

  res_ind <- get_pca_ind(pca_scaled)
  res_ind # data frame of data frames with all the info on the individuals

  ind_coord <- as.data.frame(res_ind$coord)        # Coordinates
  ind_contrib <- as.data.frame(res_ind$contrib)    # Contributions to the PCs
  ind_cos2 <- as.data.frame(res_ind$cos2)          # Quality of representation

  write_table_special(ind_coord, output_folder, file_name = "ind_coord.txt")
  write_table_special(ind_contrib, output_folder, file_name = "ind_contrib.txt")
  write_table_special(ind_cos2, output_folder, file_name = "ind_cos2.txt")

  # TODO rownames, need to change this to "variables" for output
  # TODO verify why cors are the same
  write_table_special(var_coord, output_folder, file_name = "var_coord.txt")
  write_table_special(var_cor, output_folder, file_name = "var_cor.txt")
  write_table_special(var_contrib, output_folder, file_name = "var_contrib.txt")
  write_table_special(var_cos2, output_folder, file_name = "var_cos2.txt")

  # Get pcs that are over a threshold
  bool_important <- eig_val$cumulative_variance_percent < cum_var_thres
  pcs_over_threshold <- sum(bool_important, na.rm = TRUE) + 1
  pc_name_list <- rownames(eig_val)[1:pcs_over_threshold]

  # number of items to show on scree plot
  # TODO verify what amount of items is too much to show... don't want all if
  # people are reducing all their data
  if (pcs_over_threshold >= 100) {
    pcs_to_display <- 100
  } else {
    pcs_to_display <- pcs_over_threshold
    }

  # Plot the eigenvalues (scree plot) top dimensions
  p <- fviz_eig(pca_scaled, main = "Scaled Scree Plot", choice = "variance",
         addlabels = TRUE,
         ncp = pcs_to_display)
  ggsave(p, file = paste0(output_folder, "scaled-scree-plot-top-pcs.pdf"),
         width = 14, height = 8, device = "pdf")

  # Plot the eigenvalues (scree plot) all dimensions
  p <- fviz_eig(pca_scaled, main = "Scaled Scree Plot", choice = "variance",
               addlabels = TRUE)
  ggsave(p, file = paste0(output_folder, "scaled-scree-plot-all.pdf"),
         width = 14, height = 8, device = "pdf")

  p <- fviz_pca_var(pca_scaled,
         geom_ind = "text", # ("point", "arrow", "text")
         col_var = "contrib", # color by contribs to PC
         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
         repel = TRUE
  )
  ggsave(p, file = paste0(output_folder, "variables-pca.pdf"),
         width = 14, height = 8, device = "pdf")

  # Replace original vars (remove) and add pcs
  for (pc in pc_name_list){
    pc_new_name <- paste0(pca_name, "_", pc)
    df[[pc_new_name]] <- ind_coord[[pc]]
  }
  df <- df %>% select(-all_of(feature_list))

  write_table_special(df, folder_name = output_folder,
                      file_name = "pca-dim-reduction-partial-df.txt")
  return(df)
}

all_pca_lists <- pca_list
output_folder_list <- character()
# create the folders first, so that time isn't wasted processing PCAs
for (pca_name in names(all_pca_lists)) {
  output_folder_from_group <- paste0(final_output_folder, pca_name)
  output_folder_list <- c(output_folder_list, output_folder_from_group)
}

# add main output folder to list
output_folder_list <- c(output_folder_list, final_output_folder)
# see if PCA directories exist already and if not, make them
check_for_directories(output_folder_list)

# TODO check for unique items
for (pca_name in names(all_pca_lists)) {
  df <- perform_pca(pca_list = all_pca_lists[pca_name],
                    pca_name = pca_name, df = df,
                    final_output_folder = final_output_folder)
}

# TODO write final file only after all pcs have been included loop
write_table_special(df, folder_name = final_output_folder,
                    file_name = "PCA-dim-reduction-for-workflow.txt")
