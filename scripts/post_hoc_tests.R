

#rm(list=ls()) # clear env

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


deltas_first = snakemake@input[["deltas_first"]]
deltas_previous = snakemake@input[["deltas_previous"]]
deltas_pairwise = snakemake@input[["deltas_pairwise"]]
df_raw = snakemake@input[["df_raw"]]

fixed_effects_first = snakemake@input[["fixed_effects_first"]]
fixed_effects_previous = snakemake@input[["fixed_effects_previous"]]
fixed_effects_pairwise = snakemake@input[["fixed_effects_pairwise"]]
fixed_effects_raw = snakemake@input[["fixed_effects_raw"]]

out_file = snakemake@output[["out_file"]]

response_var = snakemake@params[["response_var"]]
#absolute_values = snakemake@params[["absolute_values"]]

# notes
#exp_var = unlist(df[predictor]) res_var = unlist(df[response]) Timepoint = unlist(df["Timepoint"])

# TODO if columns dropped, allow less tests
response = response_var
#response = 'IL6'
random_effects_list = c("StudyID")
interaction = "no"

# TODO make log file containing LME errors and warnings


css_list = list(css.table = '+border:3px solid black;',
                css.caption = '+padding-top:1.2cm; padding-bottom:0.4cm',
                css.thead = '+font-weight:bold; padding:0.2cm;')

lm_post_hoc <- function(df, response, fixed_effects_list, interaction){
  fixed_effects_formula = paste(fixed_effects_list, collapse = interaction)

  reduced_formula = as.formula(paste0(response, " ~ 1"))
  full_formula = as.formula(paste0(response, " ~ ", fixed_effects_formula))

  model_reduced = lm(reduced_formula, data = df)
  model_full = lm(full_formula, data = df)
  model_anova = anova(model_full, model_reduced)

  print(model_anova)
  print(summary(model_full))
  return(model_full)
}


lme_post_hoc <- function(df, response, fixed_effects_list, random_effects_list,
                         interaction) {
  print("fixed_effects_list")
  print(fixed_effects_list)
  fixed_effects_formula = paste(fixed_effects_list, collapse = interaction)
  # linear mixed effects model (deltas)
  random_effects_formula = paste0("(1|", random_effects_list, ")",
                                   collapse = " + ")
  print("fixed_effects_formula")
  print(fixed_effects_formula)
  print("random_effects_formula")
  print(random_effects_formula)
  reduced_formula = as.formula(paste0(response, " ~ 1", " + ",
                                     random_effects_formula ))

  full_formula = as.formula(paste0(response, " ~ ", fixed_effects_formula,
                                  " + ", random_effects_formula))
  print("reduced_formula")
  print(reduced_formula)
  print("full_formula")
  print(full_formula)
  model_reduced = lmer(reduced_formula, data = df, REML=FALSE)
  model_full = lmer(full_formula, data = df, REML=FALSE)
  model_anova = anova(model_full, model_reduced)
  return(model_full)
}


if (interaction == "yes"){
  interaction = " * "
} else if (interaction == "no"){
  interaction = " + "
} else {
  stop("Interaction type incorrect")
}

# ENC_PlotName_is_P3 TODO remove modifiers

#my_tabular_model <- function()


lme_complete <- function(delta_file, fixed_effects_file, reference_type){
  html_number = 1
  reference_order_list = list("first"="01", "previous"="02", "pairwise"="03",
                              "raw" = "04")
  reference_order = reference_order_list[[reference_type]]
  df = read.csv(delta_file, sep = "\t")

  # TODO this helps find linear relns. Keep categorical?
  #df$Timepoint = factor(df$Timepoint, ordered = TRUE)

  fixed_effects_list = read.csv(fixed_effects_file, sep = "\t")
  fixed_effects_list = unique(fixed_effects_list[["decoded.features"]])

  lme_model = lme_post_hoc(df, response, fixed_effects_list,
                           random_effects_list, interaction)

  # Multivariate analysis
  tm = tab_model(lme_model,
                 file = paste0(html_number, reference_order, "-post-hoc-lme-",
                               reference_type, ".html"),
                 title = paste0("Multivariate Longitudinal Analysis - ",
                                reference_type, " as Reference"),
                 dv.labels = paste0(response, " (", reference_type,
                                    " as reference)"),
                 CSS = css_list)
  html_number = html_number + 1
  print(tm)

  lme_model2 = lme_post_hoc(df, response, c(fixed_effects_list, "Diet", "Timepoint"),
                            random_effects_list, " * ")

  # Multivariate analysis
  tm = tab_model(lme_model2,
                 file = paste0(html_number, reference_order, "-post-hoc-lme-",
                               reference_type, ".html"),
                 title = paste0("Multivariate Longitudinal Analysis ",
                                reference_type, " as Reference"),
                 dv.labels = paste0(response, " (", reference_type,
                                    " as reference)"),
                 CSS = css_list)
  html_number = html_number + 1
  print(tm)

  # Univariate analysis
  for (i in fixed_effects_list){
    lme_i = lme_post_hoc(df, response, i, random_effects_list, interaction)
    tm = tab_model(lme_i,
                   file = paste0(html_number, reference_order,
                                 "-post-hoc-lme-univariate-", reference_type,
                                 '-', response, "-by-", i, ".html"),
                   title = paste("Univariate Linear Mixed Effects Model - ",
                                 reference_type, "as reference", response, "~", i),
                   dv.labels = paste0(response, " (", reference_type,
                                      " as reference)"),
                   CSS = css_list)
    html_number = html_number + 1
    print(tm)
  }

  # Univariate with interaction *
  for (i in fixed_effects_list){
    interaction_terms = c("HIV", "Diet")
    fixed_effects_list = unique(c(i, interaction_terms))
    lme_i = lme_post_hoc(df, response, fixed_effects_list, random_effects_list, " * ")
    tm = tab_model(lme_i,
                   file = paste0(html_number, reference_order,
                                 "-post-hoc-lme-univariate-", reference_type,
                                 '-', response, "-by-", i, ".html"),
                   title = paste("Univariate Linear Mixed Effects Model w/ INTERACTION - ",
                                 reference_type, "as reference", response, "~", i),
                   dv.labels = paste0(response, " (", reference_type,
                                      " as reference)"),
                   CSS = css_list)
    html_number = html_number + 1
    print(tm)
  }
  # TODO this is returned... new function needed
  return(lme_model)
}

rt = "first"
lme_first = lme_complete(delta_file = deltas_first,
                         fixed_effects_file = fixed_effects_first,
                         reference_type = rt)

rt = "previous"
lme_previous = lme_complete(delta_file = deltas_previous,
                            fixed_effects_file = fixed_effects_previous,
                            reference_type = rt)

rt = "pairwise"
lme_pairwise = lme_complete(delta_file = deltas_pairwise,
                            fixed_effects_file = fixed_effects_pairwise,
                            reference_type = rt)


# TODO for list of models, make custom labels
######### after all RF are run
tab_model(list(lme_first, lme_previous, lme_pairwise),
          show.aic = TRUE,
          show.stat = TRUE,
          title = "Multivariate Longitudinal Analysis - All Models (note: predictors were
          not necessarily important for all models as indicated by blank cells)",
          string.stat = "Stat.",
          string.est = "Est.",
          file = "01-post-hoc-combined-first-previous-pairwise.html",
          dv.labels = c(paste(response, "(first)") ,
                        paste(response, "(previous)"),
                        paste(response, "(pairwise)")),
          CSS = css_list
)

#### IF RAW run this
rt = "raw"
lme_raw = lme_complete(delta_file = df_raw,
                       fixed_effects_file = fixed_effects_raw,
                       reference_type = rt)

tab_model(lme_raw,
          file = "9-post-hoc-raw.html",
          CSS = css_list
)

# TODO temp dir
filenames_html <- sort(Sys.glob("*.html"))
html_combine(
  combine = list(filenames_html),
  out = out_file,
)


#filenames_html <- sort(Sys.glob("*.html"))
#html_combine(
#  combine = list(filenames_html),
#  out = "test-COMBINED.html",
#)



# No longitudinal analysis, use lmes? TODO Yes, if random effects needed
#tab_model(lm_raw, title = "Cross-sectional (raw data) Analysis")
#
#tab_model(lm_raw, title = 'Cross-sectional (raw data) Analysis',
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
