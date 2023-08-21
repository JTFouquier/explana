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


add_average_shap_values <- function(df_all, ds_type){
    # Encoded features are either yes or no (1 or 0); so the average
    # shap value tells how each binary feature impacted response

    input_features <- c()
    # only binary 1 important_feature values included in plot
    enc_features <- df_all %>%
        filter(dataset == ds_type) %>%
        filter(grepl("ENC_", important_features))

    type_path <- paste0("path_", str_to_lower(ds_type))

    # all input features
    ds_path <- paste0(snakemake@config[["out"]],
    snakemake@config[[type_path]], str_to_lower(ds_type), "-input-model-df.txt")

    if (file.exists(ds_path)) {
    } else {
        return(list(df_all, input_features))
    }
    df_input <- as.data.frame(read_tsv(file = ds_path, show_col_types = FALSE))
    input_features <- colnames(df_input)

    # SHAP values dataframe
    ds_path <- paste0(snakemake@config[["out"]],
    snakemake@config[[type_path]], str_to_lower(ds_type), "-SHAP-values-df.txt")

    # continue if the model/dataset was run, or return if not
    # this will create blank columns in feature occurance plot
    if (file.exists(ds_path)) {
    } else {
        return(list(df_all, input_features))
    }

    df_shap <- as.data.frame(read_tsv(file = ds_path, show_col_types = FALSE))

    # iterate through all binary encoded features
    for (i in enc_features$important_features){

        enc_instances <- df_input %>%
            select({{i}}) %>%
            mutate(original_index = row_number()) %>%
            filter(.[[i]] == 1)

        shap_instances <- df_shap[enc_instances$original_index,]
        shap_instances <- shap_instances %>%
            select(all_of({{i}}))

        avg_shap <- round(mean(shap_instances[[i]]), digits = 2)

        df_all <- df_all %>%
            mutate(shap_val = case_when(
                important_features %in% c(i) & dataset == ds_type ~ avg_shap,
                TRUE ~ shap_val
            ))
    }

return(list(df_all, input_features))
}


make_feature_heatmaps <- function(base_path) { # nolint

    color_original <- "black"
    color_first <- "black"
    color_previous <- "black"
    color_pairwise <- "black"

    import_ds_for_summary <- function(ds_path, ds_type) {
        df <- as.data.frame(read_tsv(file = ds_path, show_col_types = FALSE))
        if ("Failed Analysis" %in% colnames(df)) {
            df <- data.frame(
            important_features = c("Failed Analysis"),
            decoded_features = c("Failed Analysis"),
            feature_importance_vals = c(10000)
            )
        }

        max_rank <- length(df$feature_importance_vals)
        df$ranked_importance <- max_rank + 1 - rank(df$feature_importance_vals)

        df[[ds_type]] <- 1
        return(df)
    }

    # TODO needs to fail gracefully when model was not run
    # as in cross-sectional studies; there would only be 'Original'
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

    df <- df_join %>%
        select(dataset, important_features) %>% # nolint
        filter(important_features != "Failed Analysis") %>%
        arrange(dataset)

    write_tsv(df,
    paste0(base_path, "all-important-features.txt"))

    # add shap values for encoded features
    df_join$shap_val <- 0

    l <- add_average_shap_values(df_all = df_join, ds_type = "Original")
    df_join <- l[[1]]
    input_features_original <- l[[2]]
    l <- add_average_shap_values(df_all = df_join, ds_type = "First")
    df_join <- l[[1]]
    input_features_first <- l[[2]]
    l <- add_average_shap_values(df_all = df_join, ds_type = "Previous")
    df_join <- l[[1]]
    input_features_previous <- l[[2]]
    l <- add_average_shap_values(df_all = df_join, ds_type = "Pairwise")
    df_join <- l[[1]]
    input_features_pairwise <- l[[2]]

    # turn into long dataframe for use with feature map
    df_join_long <- pivot_longer(df_join, cols = c("Original",
    "First", "Previous", "Pairwise"), names_to = "variable",
    values_to = "value")

    df_join_long <- df_join_long %>%
        rowwise() %>%
        mutate(figure_labels = case_when(
            value == 1 & shap_val == 0 ~ paste0(toString(ranked_importance),
            ": +\\-"),
            value == 1 ~ paste0(toString(ranked_importance), ": ",
            toString(shap_val)),
            value == 0 ~ ""
        )) %>%
        # track if feature was input but not selected; figure clarity
        mutate(is_model_input = case_when(
            value == 0 & variable == "Original"
            & important_features %in% input_features_original ~ " ",
            value == 0 & variable == "First"
            & important_features %in% input_features_first ~ " ",
            value == 0 & variable == "Previous"
            & important_features %in% input_features_previous ~ " ",
            value == 0 & variable == "Pairwise"
            & important_features %in% input_features_pairwise ~ " "
        )) %>%
        filter(important_features != "Failed Analysis") %>%
        ungroup()

    write_tsv(df_join_long,
    paste0(base_path, "df_join_long.txt"))

    # dynamically change plot height based on selected feature number
    decoded_feature_facets <- length(unique((df_join$important_features)))

    if (decoded_feature_facets < 20) {
        height_ratio <- 12 / 40
    } else {
        height_ratio <- 10.8 / 40
    }
    height <- decoded_feature_facets * height_ratio

    if (height <= 4) {
        height <- 4
    }

    ggplot(df_join_long, aes(x = factor(variable, level = c("Original",
    "First", "Previous", "Pairwise")),
    y = important_features, alpha = value,
    fill = variable, color = important_features, label = figure_labels)) +

    # labels showing included features in each model
    geom_label(label.padding = unit(0.25, "lines"),
    size = 2.2, label.size = 0.7, aes(color = is_model_input,
    label = is_model_input)) +
    scale_color_manual(values = c(" " = "#9e9c9c")) +

    # rectangular labels; white text w/ shap_values
    geom_label(color = "white", label.size = 0,
    label.padding = unit(0.25, "lines"), size = 2.9) +

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
          strip.text.y = element_text(size = 8, angle = 0),
          axis.text.x = element_text(angle = 45, hjust = 0, size = 8,
          face = "bold"),
          axis.text.y = element_text(size = 8)) +
    facet_grid(decoded_features ~ ., scales = "free_y", space = "free")

    ggsave(filename = paste0(base_path, "important-feature-occurrences.svg"),
    width = 8.5, height = height)
}

make_feature_heatmaps(base_path)