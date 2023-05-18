"value" <- "variable" <- "important_features" <- NULL
"square_border" <- "dataset" <- "horizontal_lines" <- "grid_color" <- NULL

source("scripts/install.R")
source("scripts/colors.R")
package_list <- c("ggplot2", "reshape2", "tidyverse", "readr")
install_r_packages(package_list = package_list)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(readr)


selected_features_original <- snakemake@input[["selected_features_original"]]
selected_features_first <- snakemake@input[["selected_features_first"]]
selected_features_previous <- snakemake@input[["selected_features_previous"]]
selected_features_pairwise <- snakemake@input[["selected_features_pairwise"]]

base_path <- snakemake@config[["out"]]

out_file <- snakemake@output[["out_file"]]

make_feature_heatmaps <- function(base_path) { # nolint

    color_original <- "black"
    color_first <- "black"
    color_previous <- "black"
    color_pairwise <- "black"

    import_ds_for_summary <- function(ds_path, ds_type) {
        df <- as.data.frame(read_tsv(file = ds_path))
        # if there are no rows make empty df for figure placeholder
        # need original, first, previous, pairwise
        if (dim(df)[1] == 0) {
            df <- data.frame(matrix(ncol = 3, nrow = 0))
            colnames(df) <- c("important_features", "decoded_features", ds_type)
            df$important_features <- character()
            df$decoded_features <- character()
            df[[ds_type]] <- integer()
            return(df)
        }
        df[[ds_type]] <- 1
        df <- select(df, -feature_importance_vals) # nolint
        return(df)
    }

    df_original <- import_ds_for_summary(selected_features_original, "Original")
    df_first <- import_ds_for_summary(selected_features_first, "First")
    df_previous <- import_ds_for_summary(selected_features_previous, "Previous")
    df_pairwise <- import_ds_for_summary(selected_features_pairwise, "Pairwise")

    df_original$dataset <- "Original"
    df_first$dataset <- "First"
    df_previous$dataset <- "Previous"
    df_pairwise$dataset <- "Pairwise"

    df_list <- list(df_original, df_first, df_previous, df_pairwise)
    df_join <- Reduce(function(x, y) merge(x, y, all = TRUE), df_list)

    df_join[is.na(df_join)] <- 0

    df_all <- df_join %>%
        select(dataset, important_features) # nolint

    write_tsv(df_all,
    paste0(base_path, "all-important-features.txt"))

    # dynamically change the plot height based on feature number
    decoded_feature_facets <- length(unique(df_join$important_features))
    height_per_41 <- 10
    num_features_test <- 41
    height <- decoded_feature_facets * height_per_41 / num_features_test
    ggplot(melt(df_join), aes(x = variable, y = important_features,
    alpha = value, fill = variable, color = important_features)) +
    geom_point(colour = square_border, pch = 22, size = 5) +
    scale_fill_manual(values = c("Original" = color_original,
    "First" = color_first, "Previous" = color_previous,
    "Pairwise" = color_pairwise)) +
    scale_alpha_identity(guide = "none") +
    labs(y = expression(atop("Important Features",
    atop("(Categorical variable values are encoded)"))),
    x = "Dataset") +
    scale_x_discrete(position = "top") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_line(color = grid_color),
          panel.grid.major.x = element_line(color = horizontal_lines),
          panel.grid.major.y = element_blank(),
          panel.spacing = unit(0, "lines"),
          strip.text.y = element_text(size = 11, angle = 0),
          axis.text.x = element_text(angle = 45, hjust = 0, size = 11,
          face = "bold"),
          axis.text.y = element_text(size = 11)) +
    facet_grid(decoded_features ~ ., scales = "free_y", space = "free")

    ggsave(filename = paste0(base_path, "important-feature-occurrences.svg"),
    width = 6, height = height)
}

make_feature_heatmaps(base_path)