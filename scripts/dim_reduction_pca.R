

# Script that assists with performing dimensionality reduction via
# principal component analysis (PCA) on sets of user specified input variables

# Important PCs are selected past cumulative threshold cutoff.
# The most important PCs replace the original variables in the dataframe.
#
# var1 var2 var3 might get converted to dim1 and dim2;
# Overall, a dimensionality reduction from 3 to 2 variables

# TODO verify where to clear data in project. I don't want to reload libs
# TODO check if vars exist in df, make sure no two vars are in any other list
# TODO pass height and width values
# TODO feature: use other options instead of list (slices, ranges, etc)

# rm(list=ls())
# Adapted from code written by Abigail Armstrong, modified by Jennifer Fouquier
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/#theory-behind-pca-results

library(factoextra)
library(MASS)
library(dplyr)

# user inputs from Snake
in_file = snakemake@input[["in_file"]]
out_file = snakemake@output[["out_file"]]
sample_id = snakemake@config[["sample_id"]]

pca_groups_list = snakemake@params[["pca_groups_list"]]

# TODO sort out good folder names and locations
final_output_folder = paste0(snakemake@config[["out"]],
                             snakemake@config[["path_dim_pca"]],
                             snakemake@config[["ds_name"]], "/")


write_table_special <- function(df, folder_name, file_name) {
  outf <- paste0(folder_name, file_name) # TODO
  write.table(df, out_file, row.names=FALSE, sep = "\t")
}
inf <- paste0(in_file)
df <- read.table(inf, header=TRUE, sep="\t", strip.white=FALSE)

# verify none of the PCA folders exist before doing analysis
check_for_directories <- function(folder_list){
  print(paste0("folder_list", folder_list))
  for (i in folder_list){
    if (file.exists(i)) {
      print(paste0("i", i))
      print(paste0(i, " directory already exists"))
      break
    } else {
      print(paste0("i", i))
      print(paste0(i, "created directory"))
    dir.create(i)
    }
  }
}


perform_pca <- function(pca_list, df, final_output_folder) {

  # get items from each list describing how to perform pca
  pca_name = pca_list$pca_name
  output_folder = paste0(final_output_folder, pca_name, "/")
  feature_list = pca_list$feature_list
  cum_var_thres = pca_list$cum_var_thres

  pca.df <- df %>% select(all_of(feature_list))
  # Scale PCA plot
  # TODO prcomp. this brings up need to allow user adjustments to functions
  pca.scaled <- prcomp(pca.df, scale=TRUE)

  eig.val <- get_eigenvalue(pca.scaled)
  eig.val # data frame with all info on eigen values and variance of axes

  res.var <- get_pca_var(pca.scaled)
  res.var #data frame of data frames with all the info on the variables
  var.coord <- as.data.frame(res.var$coord)     # Coordinates
  var.cor <- as.data.frame(res.var$cor)         # Corrs between vars and dims
  var.contrib <- as.data.frame(res.var$contrib) # Contributions to the PCs
  var.cos2<- as.data.frame(res.var$cos2)        # Quality of representation

  res.ind <- get_pca_ind(pca.scaled)
  res.ind #data frame of data frames with all the info on the individuals

  ind.coord <- as.data.frame(res.ind$coord)        # Coordinates
  ind.contrib <- as.data.frame(res.ind$contrib)    # Contributions to the PCs
  ind.cos2 <- as.data.frame(res.ind$cos2)          # Quality of representation

  # repetitive but not worth time
  # TODO fix these files
  write_table_special(ind.coord, output_folder, file_name="ind.coord.txt")
  write_table_special(ind.contrib, output_folder, file_name="ind.contrib.txt")
  write_table_special(ind.cos2, output_folder, file_name="ind.cos2.txt")

  # TODO rownames, need to change this to "variables" for output
  # TODO verify why cors are the same
  write_table_special(var.coord, output_folder, file_name="var.coord.txt")
  write_table_special(var.cor, output_folder, file_name="var.cor.txt")
  write_table_special(var.contrib, output_folder, file_name="var.contrib.txt")
  write_table_special(var.cos2, output_folder, file_name="var.cos2.txt")

  # TODO loop through write_table_special
  # for (i in list(var.c))

  # Get pcs that are over a threshold
  bool_important = eig.val$cumulative.variance.percent < cum_var_thres
  pcs_over_threshold = sum(bool_important, na.rm = TRUE) + 1
  pc_name_list = rownames(eig.val)[1:pcs_over_threshold]

  # number of items to show on scree plot
  # TODO verify what amount of items is too much to show... don't want all if
  # people are reducing all their data
  if (pcs_over_threshold >= 100){
    pcs_to_display = 100
  } else {
    pcs_to_display = pcs_over_threshold}

  # Plot the eigenvalues (scree plot) top dimensions
  p = fviz_eig(pca.scaled, main = "Scaled Scree Plot", choice = "variance",
         addlabels = TRUE,
         ncp=pcs_to_display)
  ggsave(p, file = paste0(output_folder, "scaled-scree-plot-top-pcs.pdf"),
         width = 14, height = 8, device='pdf')

  # Plot the eigenvalues (scree plot) all dimensions
  p = fviz_eig(pca.scaled, main = "Scaled Scree Plot", choice = "variance",
               addlabels = TRUE)
  ggsave(p, file = paste0(output_folder, "scaled-scree-plot-all.pdf"),
         width = 14, height = 8, device='pdf')

  p = fviz_pca_var(pca.scaled,
         geom.ind = "text", # type of graph ("point", "arrow", "text")
         col.var = "contrib", # color by contributions to the PC
         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # gradient colors
         repel = TRUE #avoid text overlapping
  )
  ggsave(p, file = paste0(output_folder, "variables-pca.pdf"),
         width = 14, height = 8, device='pdf')

  # Replace original vars (remove) and add pcs
  for (pc in pc_name_list){
    pc_new_name = paste0(pca_name, "_", pc)
    df[[pc_new_name]] <- ind.coord[[pc]]
  }
  df = df %>% select(-all_of(feature_list))

  # FOR REPORT

  print(paste0("pc_name_list", pc_name_list))
  print(paste0("feature_list", feature_list))

  write_table_special(df, folder_name = output_folder,
                      file_name = "pca-dim-reduction_partial_df.txt")
  return(df)
}

all_pca_lists = pca_groups_list
all_pca_lists = eval(parse(text=all_pca_lists))
output_folder_list = character()
# create the folders first, so that time isn't wasted processing PCAs
for (pca_info in all_pca_lists) {
  print(paste0("pca_info", pca_info))

  output_folder_from_group = paste0(final_output_folder, pca_info$pca_name)
  output_folder_list = c(output_folder_list, output_folder_from_group)
}

# add main output folder to list
output_folder_list = c(output_folder_list, final_output_folder)
# see if PCA directories exist already and if not, make them
check_for_directories(output_folder_list)

# TODO check for unique items
for (pca_info in all_pca_lists) {
  df = perform_pca(pca_list=pca_info, df=df,
                   final_output_folder=final_output_folder)
}

# TODO write final file only after all pcs have been included loop
write_table_special(df, folder_name = final_output_folder,
                    file_name = "final-pca-dim-reduction.txt")

