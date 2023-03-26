source("scripts/install.R")
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

    # print(df_original)
    # print(df_first)

    # df_list <- list(df_original, df_first, df_previous, df_pairwise)
    df_list <- list(df_pairwise, df_previous, df_first, df_original)
    df_join <- Reduce(function(x, y) merge(x, y, all = TRUE), df_list)

    df_join[is.na(df_join)] <- 0

    ggplot(melt(df_join), aes(x = important.features, y = variable,
    fill = decoded.features, color = decoded.features, alpha = value)) +
    geom_tile(colour = "gray50") +
    # geom_text(aes(label = )) +
    scale_alpha_identity(guide = "none") +
    coord_equal(expand = 0) +
    labs(x = expression(atop("Selected Features",
    atop(italic("Categorical variable values are encoded")))),
    fill = "Decoded\nFeatures", y = "Dataset") +
    # labs(y = "Dataset", fill = "Decoded Features",
    # x = "Selected Features (Categorical variables are encoded)") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 5.5),
          legend.title = element_text(size = 8),
          panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 55, hjust = 1, size = 6.5))

    ggsave(filename = paste0(base_path, "important-feature-occurrences.svg"))
}