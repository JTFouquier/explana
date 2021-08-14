## Variable selection in longitudinal, multiomic analytics using machine learning methods

This is an **in-progress** workflow for identification of important predictors that relate to response variables in longitudinal or cross-sectional, multiomic projects, designed specifically for microbiome research.

Omic data, especially multiomic data, is high-dimensional (more predictors than samples) which violates assumptions made by traditional statistics, often leading to inadequately used data. In order to move towards identifying causal links between microbes and a clinical response, longitudinal microbiome research is increasingly common. These longitudinal studies, where many samples are repeatedly taken from one subject over time, in an observational or interventional manner, creates additional statistical challenges because the data is not independent and clusters should be treated as a random effect. 

To address these problems, dimensionality-reduction techniques on individual omics datasets are first implemented to reduce the number of predictor variables, while at the same time reducing multicollinearity. Data from multiple sources (microbiome, immune data, clinical variables, etc.) are automatically integrated. For longitudinal datasets, a "delta dataset" is automatically created by finding all differences in feature values between timepoints, for each subject, using a specified reference point appropriate for the type of longitudinal study performed. For example, when performing a longitudinal observational study in an autism cohort to identify relationships between changes in the microbiome and changes in behavior, all intra individual pairwise comparisons for features were evaluated before performing linear mixed effects regression. Unfortunately, multiple comparison corrections in high-dimensional data possibly led to a loss of reported true relationships. In a longitudinal intervention study in HIV+ individuals, there is an expected decrease in baseline levels of systemic inflammation (IL6) after undergoing an agrarian diet intervention. This study aims to evaluate how a response variable at other time points compares to baseline response values as a reference point. Thus, different types of “delta datasets” are created as needed.

Next, predictors are ranked by their Shapley values (importance/contribution) using either a Random Forest machine learning method for fixed effects models or a Mixed Effects Random Forests (MERF) for a mix of random and fixed effects. After identifying important variables, post-hoc statistical tests (LMs or LMEs) are performed to understand, where possible, the directionality of the relationship between each predictor and the response. 

This data-driven workflow will reduce bias by allowing thorough evaluation of predictor variables and how they relate to a response variable in different study designs. While developed for longitudinal projects, through implementation of a delta dataset, this workflow can be used for cross-sectional studies. Finally, this method can be broadly used in any field through proper formatting of input files, although microbiome specific data was prioritized.

#### Software: 

I'm new to Snakemake! Please follow their tutorials. I set it up really quickly to get things going and have not explored best practices, nor is this project ready to install or clone. 

For now, I have R and Python scripts that both run when providing the name of the final pdf from the `run_random_forest` Snakemake rule which in turn uses `random_forest.py` and then the `create_deltas` Snakemake rule that uses `create_deltas.R` if needed. It steps back through the snakemake rules, finding out what it needs to build the final output file. I added the arguments to the file name for clarity using the curly braces, but this is probably not best practice. I just wanted my names to clearly say what arguments I used. 

##### Repository:

- data: original input data
- deltas: delta datasets (examples are already there, although this is an output)
- random-forest: output (examples included)
- scripts: Python and R scripts, more to come! 