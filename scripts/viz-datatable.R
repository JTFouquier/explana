# load package
rm(list=ls()) # clear env

library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(DT)
library(tidyverse)
library(dplyr)
library(colorspace)

#df = read.csv("deltas/dm-deltas-first.txt", sep = "\t")
df = read.csv("deltas/simulated-diet-deltas-pairwise.txt", sep = "\t")
out_file = "dt-viz-deltas-pairwise.html"
df.numeric = select_if(df, is.numeric)
df.other = select_if(df, negate(is.numeric))

# CATEGORICAL colors
color_palette = "ipsum"
unique_df_values = unique(as.vector(as.matrix(df.other)))
total_colors_needed = length(unique_df_values) + length(colnames(df.numeric))
multiply_colors = (total_colors_needed/length(sjplot_pal(pal = color_palette))) + 2
unique_df_colors = rep(sjplot_pal(pal = color_palette),
                       multiply_colors)[1:total_colors_needed]

# lighten to use the colors as backgrounds
unique_df_colors = lighten(
  unique_df_colors,
  amount = 0.3, # more is lighter
  method = "relative", # absolute or relative
  space = "HCL", # HCL, HLS, or combined
  fixup = TRUE
)

colors_categoric_cols = unique_df_colors[1:length(unique_df_values)]
colors_numeric_cols = unique_df_colors[length(unique_df_values):length(unique_df_colors)]

x = datatable(df, options = list(lengthMenu = c(50, 100)))

col_num = 1
for (i in colnames(df.numeric)) {
  x = x %>%
    formatStyle(i,
      background = styleColorBar(range(df.numeric[,i]), colors_numeric_cols[[col_num]]),
      backgroundSize = '98% 88%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center')

  u.col = unique(df.numeric[i])
  if (dim(u.col)[1] == 1){
    print(i)
    x = x %>%
      formatStyle(i,
      backgroundColor = styleEqual(unique(as.vector(as.matrix(df.numeric[i]))), colors_numeric_cols[[col_num]]))
  }
  col_num = col_num + 1
}

x = x %>%
  formatStyle(names(df.other),
  backgroundColor = styleEqual(unique_df_values, colors_categoric_cols))

saveWidget(x, out_file)
