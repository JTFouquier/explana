Copyright 2023 Regents of the University of Colorado. All rights reserved.

Authored by Jennifer Therese Fouquier

## Feature selection in longitudinal microbiome studies

This is an **in-progress** workflow for identification of important predictors that relate to response variables in longitudinal (and cross-sectional) microbiome studies 


High-dimensional data, where there are more predictors than samples, is increasingly common in microbiome research. High-dimensional data often violates assumptions made by traditional statistics and can lead to inadequately used data. Longitudinal studies, where many samples are repeatedly taken from one subject over time, involve additional statistical challenges because the data is not independent and clusters (often subjects) should be treated as a random effect. 


To address challenges with feature selection while working with many datasets, I am developing a computational workflow. I implement existing tools and novel methods to help generate new hypotheses via data-driven exploratory analyses. Several options for dimensionality-reduction techniques can be used on individual datasets to reduce the number of predictor variables, while also reducing multicollinearity. Data from multiple sources (microbiome, immune data, clinical variables, etc.) are automatically integrated. For longitudinal datasets, a "delta dataset" is created by finding all differences in feature values between timepoints, for each subject, using different reference points. When performing different types of longitudinal studies, such as observational or interventional, the most interesting reference point often differs. For example in interventional studies baseline values are interesting, while in an observational study, with no expected baseline response value, all pairwise comparisons might be of interest. Thus, different types of “delta datasets” are created as needed.

After obtaining delta datasets, machine learning methods, specifically Random Forest machine learning methods for fixed effects models or a Mixed Effects Random Forests (MERF) for a mix of random and fixed effects are used to find features that relate to a response. Next BorutaSHAP is implemented to rank features as important and also useful. A feature is useful if it helps the algorithm understand the response variable values better than expected by random chance or as compared to how well all shuffled features from the original dataset perform.

Next, predictors are ranked by their Shapley values (importance/contribution) with a SHAP explainer and visualized for interpretability (I.e., understanding how your variables relate to your response variable). Post-hoc statistical tests (LMs or LMEs) are then performed to understand, where possible, the directionality of the relationship between each predictor and the response, with the caveat that tree-based methods can find non-linear relationships. 

This data-driven workflow will reduce bias by allowing thorough evaluation of predictor variables and how they relate to a response variable in different longitudinal study designs. While developed for longitudinal projects, via implementation of delta datasets, this workflow can be used for cross-sectional studies. Finally, this method can be broadly used in any field through proper formatting of input files, although microbiome specific data was prioritized. 

## Software installation

#### Start with a Snakemake installation into a new Conda environment then install other software into this environment; I plan to make install easier!
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

#### Mixed effects Random Forests
https://pypi.org/project/merf/ (https://github.com/manifoldai/merf)

#### SCNIC
https://github.com/shafferm/SCNIC (with fastspar and parallel)

#### BorutaSHAP
https://pypi.org/project/BorutaShap/


## Reading
- https://towardsdatascience.com/boruta-shap-an-amazing-tool-for-feature-selection-every-data-scientist-should-know-33a5f01285c0
- https://github.com/Ekeany/Boruta-Shap
- https://shap.readthedocs.io/en/latest/
- https://christophm.github.io/interpretable-ml-book/shap.html

## Run using practice data

`snakemake --cores 6; snakemake render_report --cores 6; open HDL-Demo.html`
 