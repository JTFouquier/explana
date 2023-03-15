# in_file = snakemake@input[["in_file"]]
# out_file = snakemake@output[["out_file"]]
# pretty_print = snakemake@config[["pretty_print"]]

package_list <- c("sjPlot", "sjmisc", "sjlabelled", "DT",
"tidyverse", "dplyr", "colorspace")
package_list <- package_list[!(package_list %in%
installed.packages()[, "Package"])]
if (length(package_list) > 0) {
  install.packages(package_list)
}

library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(DT)
library(tidyverse)
library(dplyr)
library(colorspace)


build_datatable_viz <- function(df, output_file_name) {

  df <- read.csv(df, sep = "\t", check.names = FALSE)
  # TODO remove reference (change this earlier)
  df <- df[, !grepl("_reference$", names(df))]

  df_round <- function(x, digits) {
    numeric_columns <- sapply(x, mode) == "numeric"
    x[numeric_columns] <-  round(x[numeric_columns], digits)
    x
  }

  df <- df_round(df, 4)
  df.num <- select_if(df, is.numeric)
  df.other <- select_if(df, negate(is.numeric))

  # CATEGORICAL colors
  color_palette <- "ipsum"
  unique_df_vals <- unique(as.vector(as.matrix(df.other)))
  total_colors_needed <- length(unique_df_vals) + length(colnames(df.num))
  multiply_colors <- (total_colors_needed /
    length(sjplot_pal(pal = color_palette))) + 2
  unique_df_colors <- rep(sjplot_pal(pal = color_palette),
                         multiply_colors)[1:total_colors_needed]

  # lighten to use the colors as backgrounds
  unique_df_colors <- lighten(
    unique_df_colors,
    amount = 0.4, # more is lighter
    method = "relative", # absolute or relative
    space = "HCL", # HCL, HLS, or combined
    fixup = TRUE
  )

  colors_cat <- unique_df_colors[1:length(unique_df_vals)]
  colors_num <- unique_df_colors[length(unique_df_vals):
                                length(unique_df_colors)]
  # TODO this needs reorganization; negatives are not clear in the table
  x <- datatable(df, options = list(lengthMenu <- c(50, 100)))
  col_num <- 1
  for (i in colnames(df.num)) {
    x <- x %>%
      formatStyle(i,
        background = styleColorBar(range(df.num[, i]), colors_num[[col_num]]),
        # background = styleColorBar(c(0,1), colors_num[[col_num]]),
        backgroundSize = "98% 88%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center")

    u.col <- unique(df.num[i])
    if (dim(u.col)[1] == 1) {
      x <- x %>%
        formatStyle(i,
        backgroundColor = styleEqual(unique(as.vector(
          as.matrix(df.num[i]))), colors_num[[col_num]]))
    }
    col_num <- col_num + 1
  }
  x <- x %>%
    formatStyle(names(df.other),
    backgroundColor = styleEqual(unique_df_vals, colors_cat))
  saveWidget(x, output_file_name)
}

# build_datatable(in_file, output_file_name=out_file)