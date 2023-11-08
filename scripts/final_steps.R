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


add_average_shap_values <- function(df, ds) {
    # Encoded features are either yes or no (1 or 0); so the average
    # shap value tells how each binary feature impacted response

    input_features <- c()
    # only binary 1 important_feature values included in plot
    enc_features <- df %>%
        filter(dataset == ds) %>%
        filter(grepl("ENC_", important_features))

    type_path <- paste0("path_", stringr::str_to_lower(ds))

    # all input features
    path_prefix <- paste0(snakemake@config[["out"]],
    snakemake@config[[type_path]], stringr::str_to_lower(ds))

    fp <- paste0(path_prefix, "-input-model-df.txt")
    if (file.exists(fp)) {
    } else {
        return(list(df, input_features))
    }
    df_input <- read_tsv(file = fp, show_col_types = FALSE)
    input_features <- input_features <- colnames(df_input)

    fp <- paste0(path_prefix, "-SHAP-values-df.txt")
    # continue if the model/dataset was run, or return if not
    # this will create blank columns in feature occurance plot
    if (file.exists(fp)) {
    } else {
        return(list(df, input_features))
    }
    df_shap <- as.data.frame(read_tsv(file = fp, show_col_types = FALSE))

    # iterate through all binary encoded features
    for (i in enc_features$important_features){

        enc_instances <- df_input %>%
            dplyr::select({{i}}) %>%
            dplyr::mutate(original_index = dplyr::row_number()) %>%
            filter(.[[i]] == 1)

        shap_instances <- df_shap[enc_instances$original_index, ]
        shap_instances <- shap_instances %>%
            dplyr::select(dplyr::all_of({{i}}))

        avg_shap <- round(mean(shap_instances[[i]]), digits = 2)

        df <- df %>%
            dplyr::mutate(shap_val = dplyr::case_when(
                important_features %in% c(i) & dataset == ds ~ avg_shap,
                TRUE ~ shap_val
            ))
    }

return(list(df, input_features))
}


make_summary_log_dataframe <- function() {
    load_log_df <- function(ds) {
        type_path <- paste0("path_", stringr::str_to_lower(ds))
        path_prefix <- paste0(snakemake@config[["out"]],
        snakemake@config[[type_path]], stringr::str_to_lower(ds))

        fp <- paste0(path_prefix, "-log-df.txt")
        df <- as.data.frame(read_tsv(file = fp))
        return(df)
    }

    df_original <- load_log_df("Original")
    df_first <- load_log_df("First")
    df_previous <- load_log_df("Previous")
    df_pairwise <- load_log_df("Pairwise")

    df_list <- list(df_original, df_first, df_previous, df_pairwise)
    df <- Reduce(function(x, y) merge(x, y, all = TRUE), df_list)

    write_tsv(df, paste0(base_path, "summary-log-table.txt"))

}


main <- function(base_path) { # nolint

    color_original <- "black"
    color_first <- "black"
    color_previous <- "black"
    color_pairwise <- "black"

    make_summary_log_dataframe()

    import_ds_for_summary <- function(ds_path, ds) {
        df <- as.data.frame(read_tsv(file = ds_path, show_col_types = FALSE))
        if ("Failed Analysis" %in% colnames(df)) {
            df <- data.frame(
            important_features = c("Failed Analysis"),
            decoded_features = c("Failed Analysis"),
            feature_importance_vals = c(10000)
            )
        }

        max_rank <- length(df$feature_importance_vals)
        df$ranked_importance <- max_rank + 1 - rank(df$feature_importance_vals,
        ties.method = "first")

        df[[ds]] <- 1
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

    df <- df_join %>%
        select(dataset, important_features) %>% # nolint
        filter(important_features != "Failed Analysis") %>%
        arrange(dataset)

    write_tsv(df,
    paste0(base_path, "all-important-features.txt"))

    # add shap values for encoded features
    df_join$shap_val <- 0

    l <- add_average_shap_values(df = df_join, ds = "Original")
    df_join <- l[[1]]
    input_features_original <- l[[2]]
    l <- add_average_shap_values(df = df_join, ds = "First")
    df_join <- l[[1]]
    input_features_first <- l[[2]]
    l <- add_average_shap_values(df = df_join, ds = "Previous")
    df_join <- l[[1]]
    input_features_previous <- l[[2]]
    l <- add_average_shap_values(df = df_join, ds = "Pairwise")
    df_join <- l[[1]]
    input_features_pairwise <- l[[2]]

    # turn into long dataframe for use with feature map
    df_join_long <- tidyr::pivot_longer(df_join, cols = c("Original",
    "First", "Previous", "Pairwise"), names_to = "variable",
    values_to = "value")

    print(df_join_long)

    df_join_long <- df_join_long %>%
        dplyr::rowwise() %>%
        dplyr::mutate(figure_labels = dplyr::case_when(value == 1 &
        shap_val == 0 ~ paste0(toString(ranked_importance),
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
        filter(important_features != "no_selected_features") %>%
        dplyr::ungroup()

    # get the total number of features selected by dataset for figure
    original_n <- df_join_long %>% filter(dataset == "Original")
    first_n <- df_join_long %>% filter(dataset == "First")
    previous_n <- df_join_long %>% filter(dataset == "Previous")
    pairwise_n <- df_join_long %>% filter(dataset == "Pairwise")

    all_selected_features <- df_join_long %>%
        dplyr::select(important_features, variable, value) %>%
        filter(value == 1)


    write_tsv(all_selected_features, paste0(base_path,
    "selected_features_for_fscore.txt"))

    # all input features for f-score analysis for testing
    all_input_features <- unique(c(input_features_original,
    input_features_first, input_features_previous, input_features_pairwise))
    input_feature_df <- data.frame(input_features = all_input_features)
    write_tsv(input_feature_df, paste0(base_path,
    "input_features_for_fscore.txt"))

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

    # long var names smoosh graph (esp. for taxonomy); truncate if needed
    trunc_label_len <- 40  # must be even number

    # do this for feature names
    df_join_long$temp_label <- ifelse(nchar(df_join_long$important_features) <
    trunc_label_len, df_join_long$important_features,
    paste(substr(df_join_long$important_features, 1, trunc_label_len / 2),
    "...", substr(df_join_long$important_features,
    nchar(df_join_long$important_features) - ((trunc_label_len / 2) - 1),
    nchar(df_join_long$important_features))))

    # and for decoded feature names
    df_join_long$temp_label_decoded <- ifelse(nchar(df_join_long$decoded_features) < trunc_label_len,
    df_join_long$decoded_features,
    paste(substr(df_join_long$decoded_features, 1, trunc_label_len / 2),
    "...", substr(df_join_long$decoded_features,
    nchar(df_join_long$decoded_features) - ((trunc_label_len / 2) - 1),
    nchar(df_join_long$decoded_features))))

    get_n <- function(df) {
        n_string <- toString(max(df$ranked_importance))
        if (n_string == "-Inf") {
            n_string <- ""
        } else {
            n_string <- paste0("\n   N:", n_string)
        }
        return(n_string)
    }

    custom_labels <- c(
    "Original" = paste0("Original", get_n(original_n)),
    "First" = paste0("First", get_n(first_n)),
    "Previous" = paste0("Previous", get_n(previous_n)),
    "Pairwise" = paste0("Pairwise", get_n(pairwise_n))
    )

    ggplot(df_join_long, aes(x = factor(variable, level = c("Original",
    "First", "Previous", "Pairwise")), y = temp_label, alpha = value,
    fill = variable, color = temp_label, label = figure_labels)) +

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
    scale_x_discrete(position = "top", labels = custom_labels) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_line(color = grid_color),
          panel.grid.major.x = element_line(color = horizontal_lines),
          panel.grid.major.y = element_blank(),
          panel.spacing = unit(0, "lines"),
          strip.text.y = element_text(size = 7, angle = 0),
          axis.text.x = element_text(angle = 45, hjust = 0, size = 7,
          face = "bold"),
          axis.text.y = element_text(size = 7)) +
    facet_grid(temp_label_decoded ~ ., scales = "free_y", space = "free")

    ggsave(filename = paste0(base_path, "important-feature-occurrences.svg"),
    width = 6.5, height = height)

    # remove temporary truncated label names (for fig)
    df_join_long$temp_label <- NULL
    df_join_long$temp_label_decoded <- NULL

    write_tsv(df_join_long,
    paste0(base_path, "df_join_long.txt"))

}

main(base_path)