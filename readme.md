**Copyright 2023-2024 Regents of the University of Colorado. All rights reserved.**

Authored by Jennifer Therese Fouquier

**Free for academic use.** For other licensing information and user guide can be found at [explana.io](https://www.explana.io)

Please contact me with any issues, suggestions or questions at jennifer@explana.io

## EXPLANA: A user-friendly workflow for EXPLoratory ANAlysis and feature selection in longitudinal and cross-sectional microbiome studies

EXPLANA was developed to streamline identification of important features that relate to response variables in longitudinal, and cross-sectional, microbiome studies.

Response variables can be numerical or binary categorical and input features can be numerical or categorical. If order of categorical value changes impacts a response, this will be identified.

Features are ranked by their importance, and feature impact on response is additinoally provided. An interactive .html report is generated for each analysis that visually and textually summarizes results.

**This is an exploratory method and this should be made clear when communicating results.**

# Install

### Source code

Run the command:

> `git clone git@github.com:JTFouquier/explana.git`

Next, navigate to directory with `Snakefile` and work from there.

- - -

### Create Conda environment using Mamba

You must be inside the directory with `Snakefile` to create the `explana` environment and run your analyses.

You will need to install Mambaforge [(download Mambaforge here)](https://github.com/conda-forge/miniforge#mambaforge)

> `mamba env create -f conda_envs/environment.yaml -n explana`

Activate the `explana` conda environment. This is needed to ensure specific software versions are used.

> `conda activate explana`

- - -

### Test install

Run EXPLANA with small install demo. If possible, use more cores (I use 6).

> `snakemake --configfile config/config-install-test.yaml --cores 2 --use-conda`

Open html report (or double click it from within directory)

> `open workflow-results/EXPLANA-install-test/EXPLANA-report.html`

To retest, remove output folder.

To remove and retest (`rm -r` removes a folder and *cannot* be undone):
> `rm -r workflow-results/EXPLANA-install-test/`

If workflow partially runs, rerun `snakemake` command and Snakemake will resume analysis.

- - -

## FAQs and user manual coming soon

## Useful Reading

* https://towardsdatascience.com/boruta-shap-an-amazing-tool-for-feature-selection-every-data-scientist-should-know-33a5f01285c0
* https://github.com/Ekeany/Boruta-Shap
* https://shap.readthedocs.io/en/latest/
* https://christophm.github.io/interpretable-ml-book/shap.html