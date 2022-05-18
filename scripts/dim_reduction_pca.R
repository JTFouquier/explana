

# Script that assists with performing dimensionality reduction via
# principal component analysis (PCA) on sets of user input variables

# Important PCs are selected past cumulative threshold cutoff.
# The most important PCs replace the original variables in the dataframe.
#
# var1 var2 var3 might get converted to dim1 and dim2; reduction of variables
# from 3 to 2.

rm(list=ls())
# Adapted from code written by Abigail Armstrong, modified by Jennifer Fouquier
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/#theory-behind-pca-results

library(factoextra)
library(MASS)
library(dplyr)

setwd("/Users/jenniferfouquier/repos/snakemake-merf/")

# TODO will later be user inputs
variance_threshold = 80
input_folder = "data/"
input_file = "real-data-no-asvs.txt"
output_folder = "data/pca_output/"

feature_list1 <- c('Triglycerides', 'LDL', 'Leptin', 'Adiponectin',
                   'Insulin', 'Glucose')

# use this to rename pcs for clarity
# name_pcs = c("pca1")
# group_1 = c("pca1", feature_list1, 70)
# group_2 = c("pca2", feature_list2, 90)

write_table_special <- function(df, folder_name, file_name) {
  outf <- paste(folder_name, file_name, sep ="")
  df = data.frame("StudyID.Timepoint"=rownames(df), df)
  write.matrix(as.matrix(df), file = outf, sep='\t')
}

inf <- paste(input_folder, input_file, sep="")
df <- read.table(inf, header=TRUE, sep="\t", strip.white=FALSE,
                 row.names='SampleID')
pca.df <- subset(df, select=feature_list1)


# Scale PCA plot
# TODO prcomp. this brings up need to allow user adjustments to functions
pca.scaled <- prcomp(pca.df, scale=TRUE)

eig.val <- get_eigenvalue(pca.scaled)
eig.val # data frame with all info on eigen values and variance of axes

res.axis.perc <- eig.val$variance.percent # get percent variation explained for each axis

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

# TODO this needs to be fixed. loop
# dfs_to_write <- list(ind.coord, ind.contrib, ind.cos2)
# df_file_names <- c("ind.coord.txt", "ind.contrib.txt", "ind.cos2.txt")
# for(i in seq_along(dfs_to_write)){
#   print(i)
#   print(typeof(dfs_to_write[i][0]))
#   write_table_special(df=dfs_to_write[i], folder_name=output_folder,
#                       file_name=df_file_names[i])
# }

# repetitive, use code above after fixing
write_table_special(ind.coord, output_folder, file_name="ind.coord.txt")
write_table_special(ind.contrib, output_folder, file_name="ind.contrib.txt")
write_table_special(ind.cos2, output_folder, file_name="ind.cos2.txt")

# TODO rownames, need to change this to "variables" for output
write_table_special(var.coord, output_folder, file_name="var.coord.txt")
write_table_special(var.cor, output_folder, file_name="var.cor.txt")
write_table_special(var.contrib, output_folder, file_name="var.contrib.txt")
write_table_special(var.cos2, output_folder, file_name="var.cos2.txt")

# Get pcs that are over a threshold
bool_important = eig.val$cumulative.variance.percent < variance_threshold
pcs_over_threshold = sum(bool_important, na.rm = TRUE) + 1
pc_name_list = rownames(eig.val)[1:pcs_over_threshold]

# number of items to show on scree plot
if (pcs_over_threshold >= 100){
  pcs_to_display = 100
} else {
  pcs_to_display = pcs_over_threshold}


# Plot the eigenvalues (scree plot)
pdf(file=paste(output_folder, "scaled-scree-plot-top-pcs.pdf", sep=""),
    height=8, width=14)
fviz_eig(pca.scaled, main = "Scaled Scree Plot", choice = "variance",
         addlabels = TRUE,
         ncp=pcs_to_display)
dev.off()

pdf(file=paste(output_folder, "scaled-scree-plot-all.pdf", sep=""), height=8,
    width=14)
fviz_eig(pca.scaled, main = "Scaled Scree Plot", choice = "variance",
         addlabels = TRUE)
dev.off()

pdf(file=paste(output_folder, "variables-pca.pdf", sep=""), height=8, width=14)
fviz_pca_var(pca.scaled,
             geom.ind = "text", #specifies the type of graph (i.e. "point", "arrow", or "text")
             col.var = "contrib", #color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), #specifying the gradient colors
             repel = TRUE #avoid text overlapping
)
dev.off()


# Replace original vars (remove) and add pcs
df_updated = df
for (pc in pc_name_list){
  df_updated[[pc]] <- ind.coord[[pc]]
}
df_updated = df_updated %>% dplyr::select(!feature_list1)

# TODO write file only after all pcs have been included loop
write_table_special(df_updated, folder_name = output_folder,
                    file_name = "updated_df_dim_reduction.txt")

