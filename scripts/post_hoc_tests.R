source("scripts/install.R")
package_list <- c("dplyr", "R3port", "usedist", "lmerTest", "lme4",
"sjPlot", "sjmisc", "sjlabelled", "ggpubr", "ggplot2", "stringr", "pdftools")
install_r_packages(package_list = package_list)

library(dplyr)
library(R3port)
library(usedist)
library(lmerTest)
library(lme4)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(ggpubr)
library(ggplot2)
library(stringr)
library(pdftools)

source("scripts/post_hoc_visualizations.R")


# TODO use dataset after drop for RF vs before
deltas_first <- snakemake@input[["deltas_first"]]
deltas_previous <- snakemake@input[["deltas_previous"]]
deltas_pairwise <- snakemake@input[["deltas_pairwise"]]
df_original <- snakemake@input[["df_original"]]

important_features_first <- snakemake@input[["fixed_effects_first"]]
important_features_previous <- snakemake@input[["fixed_effects_previous"]]
important_features_pairwise <- snakemake@input[["fixed_effects_pairwise"]]
important_features_original <- snakemake@input[["fixed_effects_original"]]

out_file <- snakemake@output[["out_file"]]

output_folder <- paste0(snakemake@config[["out"]],
                       snakemake@config[["path_post_hoc"]])

response <- snakemake@config[["response_var"]]
random_effect <- snakemake@config[["random_effect"]]
timepoint <- snakemake@config[["timepoint"]]

#absolute_values <- snakemake@params[["absolute_values"]]

# notes
#exp_var <- unlist(df[predictor]) res_var = unlist(df[response]) Timepoint = unlist(df["Timepoint"])

# TODO if columns dropped, allow less tests
# TODO get control/references vars for interactions
interaction <- "no"
# TODO make log file containing LME errors and warnings

dir.create(output_folder)

my_css <- function(analysis_type, dataset_type) {

  # Outside border of table
  border_color_list <- list("multivariate" = "+border: 18px solid #CAE1FF; margin-left: auto; margin-right: auto;",
                            "univariate" = "+border: 18px solid #FFE1FF; margin-left: auto; margin-right: auto;")
  border_css <- border_color_list[[analysis_type]]

  # Title banner color
  banner_color <- function(banner_color) {
    i <- paste("+padding: 14px; background-color:", banner_color,
               "; border-top: 36px solid white; font-style:normal; font-size: 18px; ")
    return(i)
  }
  banner_color_list <- list("first" = banner_color("cornflowerblue"),
                            "previous" = banner_color("lightgreen"),
                            "pairwise" = banner_color("#AB82FF"),
                            "original" = banner_color("#FFFACD"),
                            "combined" = banner_color("#BABABA"))
  banner_css <- banner_color_list[[dataset_type]]

  css_list <- list(css.table = border_css,
                   css.body = "+border: 10px solid; margin-left: auto;
                   margin-right: auto; padding: 12px;",
                   css.caption = banner_css,
                   css.thead = "+font-weight:bold; font-style:normal;
                   font-size: 15px; padding-top: 9px; padding-bottom: 9px; 
                   background-color: #F2F2F2;",
                   css.summary = "+background-color: #F2F2F2;",
                   css.depvarhead = "+color: black; padding:5px;",
                   css.randomparts = "+font-weight:bold; text-align:left;
                   padding-top:.8em; background-color: #F2F2F2;")
  return(css_list)
}

lm_post_hoc <- function(df, response, fixed_effects_list, interaction) {
  fixed_effects_formula <- paste(fixed_effects_list, collapse = interaction)

  reduced_formula <- as.formula(paste0(response, " ~ 1"))
  full_formula <- as.formula(paste0(response, " ~ ", fixed_effects_formula))
  model_reduced <- lm(reduced_formula, data = df)
  model_full <- lm(full_formula, data = df)
  model_anova <- anova(model_full, model_reduced)

  print(model_anova)
  print(summary(model_full))
  return(model_full)
}


lme_post_hoc <- function(df, response, important_features_list, random_effect,
                         interaction) {

  important_features_formula <- paste0("`", important_features_list, "`",
                                      collapse = interaction)
  random_effects_formula <- paste0("(1|", random_effect, ")",
                                  collapse = " + ")

  reduced_formula <- as.formula(paste0(response, " ~ 1", " + ",
                                      random_effects_formula))

  full_formula <- as.formula(paste0(response, " ~ ", important_features_formula,
                                   " + ", random_effects_formula))

  model_reduced <- lmer(reduced_formula, data = df, REML = FALSE)
  model_full <- lmer(full_formula, data = df, REML = FALSE)
  model_anova <- anova(model_full, model_reduced)
  return(model_full)
}


if (interaction == "yes") {
  interaction <- " * "
} else if (interaction == "no") {
  interaction <- " + "
} else {
  stop("Interaction type incorrect")
}


lme_complete <- function(df_file, important_features_file, reference_type) {
  html_number <- 100
  reference_order_list <- list("original" = "01", "first" = "02",
                               "previous" = "03", "pairwise" = "04")
  reference_order <- reference_order_list[[reference_type]]
  df <- read.csv(df_file, sep = "\t", check.names = FALSE)
  # TODO this helps find linear relns. Keep categorical?
  df[[paste0(timepoint)]] <- factor(df[[paste0(timepoint)]], ordered = TRUE)

  important_features_list <-
  read.csv(important_features_file, sep = "\t", check.names = FALSE)
  important_features_list <-
  unique(important_features_list[["decoded.features"]])

  # Multivariate analysis (all important variables explain response)
  css_list <- my_css(analysis_type = "multivariate",
  dataset_type <- reference_type)
  lme_model <- lme_post_hoc(df, response, important_features_list,
                           random_effect, interaction)
  print(tab_model(lme_model,
                  file = paste0(output_folder, reference_order,
                                "-", html_number,
                                "-multivariate-individual-",
                                reference_type, ".html"),
                  title = paste0("Multivariate Model - ", reference_type),
                  dv.labels = paste0(response, " (", reference_type, ")"),
                  CSS = css_list))
  html_number <- html_number + 1

  # Univariate analysis
  css_list <- my_css(analysis_type = "univariate",
  dataset_type = reference_type)
  for (i in important_features_list){
    lme_i <- lme_post_hoc(df, response, i, random_effect, interaction)
    print(tab_model(lme_i,
                    file = paste0(output_folder, reference_order,
                                  "-", html_number,
                                  "-univariate-", reference_type,
                                  "-", response, "-by-", i, ".html"),
                    title = paste("Univariate Model - ", reference_type,
                                  response, "~", i),
                    dv.labels = paste0(response, " (", reference_type, ")"),
                    CSS = css_list))
    html_number <- html_number + 1
  }
  return(lme_model)
}


# Do all statistics for each type of dataset
lme_first <-
lme_complete(df_file = deltas_first,
             important_features_file = important_features_first,
             reference_type = "first")
lme_previous <-
lme_complete(df_file = deltas_previous,
             important_features_file = important_features_previous,
             reference_type = "previous")
lme_pairwise <-
lme_complete(df_file = deltas_pairwise,
             important_features_file = important_features_pairwise,
             reference_type = "pairwise")
lme_original <-
lme_complete(df_file = df_original,
             important_features_file = important_features_original,
             reference_type = "original")

# Combine delta results for simplified interpretation
tab_model(list(lme_first, lme_previous, lme_pairwise),
          show.aic = TRUE,
          show.stat = TRUE,
          title = "Multivariate Analysis - Delta Dataset Models",
          string.stat = "Stat.",
          string.est = "Est.",
          file = paste0(output_folder, "01-multivariate-deltas.html"),
          dv.labels = c(paste(response, "(first)"),
                        paste(response, "(previous)"),
                        paste(response, "(pairwise)")),
          CSS = my_css("multivariate", "combined")
)

filenames_html_multivariate_combined <-
sort(Sys.glob(file.path(output_folder, "*multivariate*deltas*.html")))
filenames_html_multivariate_individual <-
sort(Sys.glob(file.path(output_folder, "*multivariate*individual*.html")))
filenames_html_univariate <-
sort(Sys.glob(file.path(output_folder, "*univariate*.html")))

filenames_html_complete <- c(filenames_html_multivariate_combined,
                            filenames_html_multivariate_individual,
                            filenames_html_univariate)

html_combine(
  combine = list(filenames_html_complete),
  out = out_file,
)

update_pdf_list <- function(pdf_list, file_to_add) {
  if (file.exists(file_to_add)) {
    pdf_list <- c(pdf_list, file_to_add)
  }
 return(pdf_list)
}


pdf_list <- c()

# dataframes for all datasets to run post hoc tests on
df_list <- c(df_original, deltas_first, deltas_previous, deltas_pairwise)

# top features selected from each dataset
important_features_list <- c(important_features_original,
important_features_first, important_features_previous,
important_features_pairwise)

# paths for output for each dataset
df_path_list <- c(snakemake@config[["path_rf_original"]],
snakemake@config[["path_rf_first"]], snakemake@config[["path_rf_previous"]],
snakemake@config[["path_rf_pairwise"]])

df_type_list <- c("Original", "First", "Previous", "Pairwise")

for (i in seq_along(df_list)) {
  out_file <- post_hoc_viz(df = df_list[i],
  important_features = important_features_list[i],
  gg_info = df_type_list[i], df_path = df_path_list[i])

  # if the file was created, add it to list of
  # pdfs for original, first, previous, pairwise
  # datasets
  pdf_list <- update_pdf_list(pdf_list, out_file)
}

pdf_combine(pdf_list, output = paste0(snakemake@config[["out"]],
"post-hoc-combined.pdf"))



# No longitudinal analysis, use lmes? TODO Yes, if random effects needed
#tab_model(lm_original, title = "Cross-sectional (original data) Analysis")
#
#tab_model(lm_original, title = 'Cross-sectional (original data) Analysis',
#           CSS = list(css.table = 'border:2px solid red;',
#                     css.summary = 'font-weight:bold;',
#                     css.lasttablerow = 'border-bottom: 1px dotted blue;',
#                     css.colnames = 'color:green;',
#                     css.arc = 'color:blue;',
#                     css.caption = 'color:blue;')
#)

# TODO model (raw is rank deficient)
# TODO Previous model is rank deficient!
# TODO ? Random effects results "StudyID

# legacy code and other TODOs:

# how to know which one is the reference level when I don't select it
#df$DailyWeatherOfficial <- factor(df$DailyWeatherOfficial,
# levels = c("Windy_Windy", "Microburst_Windy", "Windy_Microburst"))

# TODO add beta diversity merge to metadata for easier access later
# TODO df is what file... merge data before
# TODO rescale vars
# verify that order works when 1_2, 1_4, 1_3 for characters... ? sort; levels
# https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_model_estimates.html

### TODO figure out if references and original values are necessary
#reference_is_Windy
#DailyWeatherOfficial_reference
#ResponseMint_reference
#FertilizerLBrand_reference_is_B2
#FertilizerLBrand (was not included, just reference)

# TODO ideas for making viz better
# TODO highlight rare things using CSS for deltas and raw values
# TODO add dots and borders in between StudyIDs?
# figure out CSS for displaying interesting things in dfs
# gradients? categorical, etc.

# Example lme
# model_reduced = lmer(HDL ~ 1 + (1|StudyID), data = df, REML=FALSE)
# model_full = lmer(HDL ~ Timepoint + (1|StudyID), data = df, REML=FALSE)
# print(summary(model_full))