source("scripts/install.R")
source("scripts/colors.R")
package_list <- c("ggplot2", "reshape2")
install_r_packages(package_list = package_list)
library(ggplot2)
library(reshape2)


make_feature_heatmaps <- function(base_path) {

    import_dataset_for_summary <- function(ds_type) {
        df <- as.data.frame(read_tsv(file = paste0(base_path,
            "04-SELECTED-FEATURES-", ds_type, "/", ds_type,
            "-boruta-important.txt")))
        # if there are no rows make empty df for figure placeholder
        # need original, first, previous, pairwise
        if (dim(df)[1] == 0) {
            df <- data.frame(matrix(ncol = 3, nrow = 0))
            colnames(df) <- c("important.features", "decoded.features", ds_type)
            df$important.features <- character()
            df$decoded.features <- character()
            df[[ds_type]] <- integer()
            return(df)
        }
        df[[ds_type]] <- 1
        df <- select(df, -feature_importance_vals)
        return(df)
    }

    df_original <- import_dataset_for_summary("original")
    df_first <- import_dataset_for_summary("first")
    df_previous <- import_dataset_for_summary("previous")
    df_pairwise <- import_dataset_for_summary("pairwise")

    df_list <- list(df_original, df_first, df_previous, df_pairwise)
    df_join <- Reduce(function(x, y) merge(x, y, all = TRUE), df_list)

    df_join[is.na(df_join)] <- 0

    # dynamically change the plot height based on feature number
    decoded_feature_facets <- length(unique(df_join$important.features))
    height_per_41 <- 9.4
    num_features_test <- 41
    height <- decoded_feature_facets * height_per_41 / num_features_test

    ggplot(melt(df_join), aes(x = variable, y = important.features,
    alpha = value, fill = variable, color = important.features)) +
    geom_point(colour = square_border, pch = 22, size = 5.8) +
    scale_fill_manual(values = c("original" = color_original,
    "first" = color_first, "previous" = color_previous,
    "pairwise" = color_pairwise)) +
    scale_alpha_identity(guide = "none") +
    labs(y = expression(atop("Important Features",
    atop(italic("Categorical variable values are encoded")))),
    x = "Dataset") +
    scale_x_discrete(position = "top", labels = c("original" = "Original",
    "first" = "First", "previous" = "Previous", "pairwise" = "Pairwise")) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_line(color = grid_color),
          panel.grid.major.x = element_line(color = horizontal_lines),
          panel.grid.major.y = element_blank(),
          panel.spacing = unit(0, "lines"),
          strip.text.y = element_text(size = 8, angle = 0),
          axis.text.x = element_text(angle = 45, hjust = 0, size = 8,
          face = "bold"),
          axis.text.y = element_text(size = 8)) +
    facet_grid(decoded.features ~ ., scales = "free_y", space = "free")

    ggsave(filename = paste0(base_path, "important-feature-occurrences.svg"),
    width = 6, height = height)
}