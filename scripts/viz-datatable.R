# load package

library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(DT)
library(tidyverse)
library(dplyr)
library(colorspace)

df = read.csv("deltas/dm-deltas-first.txt", sep = "\t")
df = read.csv("deltas/simulated-plants-deltas-previous.txt", sep = "\t")

df$Timepoint = factor(df$Timepoint, ordered = TRUE)

df.numeric = select_if(df, is.numeric)

# TODO this is hardcoded, needs to be fixed for different datasets
df.numeric0 <- df.numeric %>% select_if(~max(unique(as.vector(as.matrix(.))))<=1)
df.numeric1 <- df.numeric %>% select_if(~max(unique(as.vector(as.matrix(.))))>1 && max(unique(as.vector(as.matrix(.))))<=10)
df.numeric2 <- df.numeric %>% select_if(~max(unique(as.vector(as.matrix(.))))>10 && max(unique(as.vector(as.matrix(.))))<=100)
df.numeric3 <- df.numeric %>% select_if(~max(unique(as.vector(as.matrix(.))))>100 && max(unique(as.vector(as.matrix(.))))<=1000)
df.numeric4 <- df.numeric %>% select_if(~max(unique(as.vector(as.matrix(.))))>1000 && max(unique(as.vector(as.matrix(.))))<=10000)

names(df.numeric1)
names(df.numeric0)
names(df.numeric1)
names(df.numeric2)
names(df.numeric3)
names(df.numeric)

df.other = select_if(df, negate(is.numeric))

color_palette = "ipsum"
unique_df_values = unique(as.vector(as.matrix(df.other)))
multiply_colors = length(unique_df_values)/length(sjplot_pal(pal = color_palette)) + 2
unique_df_colors = rep(sjplot_pal(pal = color_palette),
                       multiply_colors)[1:length(unique_df_values)]

# lighten to use the colors as backgrounds
unique_df_colors = lighten(
  unique_df_colors,
  amount = 0.3, # more is lighter
  method = "relative", # absolute or relative
  space = "HCL", # HCL, HLS, or combined
  fixup = TRUE
)

x = datatable(df) %>%
  formatStyle(names(df.numeric0),
  background = styleColorBar(range(df.numeric0), unique_df_colors[[6]] ),
  backgroundSize = '98% 88%',
  backgroundRepeat = 'no-repeat',
  backgroundPosition = 'center') %>%

  formatStyle(names(df.numeric1),
  background = styleColorBar(range(df.numeric1), color = unique_df_colors[[1]]),
  backgroundSize = '98% 88%',
  backgroundRepeat = 'no-repeat',
  backgroundPosition = 'center') %>%

  formatStyle(names(df.numeric2),
  background = styleColorBar(range(df.numeric2), color = unique_df_colors[[2]]),
  backgroundSize = '98% 88%',
  backgroundRepeat = 'no-repeat',
  backgroundPosition = 'center') %>%

  formatStyle(names(df.numeric3),
  background = styleColorBar(range(df.numeric3), color = unique_df_colors[[3]]),
  backgroundSize = '98% 88%',
  backgroundRepeat = 'no-repeat',
  backgroundPosition = 'center') %>%

  formatStyle(names(df.numeric4),
  background = styleColorBar(range(df.numeric4), color = unique_df_colors[[4]]),
  backgroundSize = '98% 88%',
  backgroundRepeat = 'no-repeat',
  backgroundPosition = 'center') %>%


  formatStyle(names(df.other),
  backgroundColor = styleEqual(unique_df_values, unique_df_colors))

saveWidget(x, 'dt-table.html')