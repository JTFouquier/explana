**Copyright 2023-2024 Regents of the University of Colorado. All rights reserved.**

Authored by Jennifer Therese Fouquier

EXPLANA: A user-friendly workflow for EXPLoratory ANAlysis and feature selection in longitudinal and cross-sectional microbiome studies
========================================================================

.. image:: https://github.com/JTFouquier/explana/blob/main/images/report-screenshot.png

EXPLANA was developed to streamline identification of important features that relate to outcome/response variables in longitudinal, and cross-sectional, microbiome studies.

Outcome variables can be numerical or binary categorical and input features can be numerical or categorical. If order of categorical value changes impact a response, this will be identified.

Features are ranked by their importance, and feature impact on response is provided. An interactive .html report is generated for each analysis that visually and textually summarizes results.

**This is an exploratory method and this should be made clear when communicating results.**

User Guide
===========

`explana.io/documentation<https://www.explana.io/documentation>`_

Using EXPLANA
==============

You can use this tool in several ways; it will be available under multiple licenses:

1) **Free for academic use.** Academic users can install and use the software.
2) **Analytic service** I am happy to analyze your data for you as a paid service. We would discuss your project, I would analyze your data, then follow up with a meeting to go over the reports/results. For more information, see `www.jenniferfouquier.com<https://www.jenniferfouquier.com>`_ or `schedule time to chat<https://www.jenniferfouquier.com/appointments>`_. No charge for consultation.
3) **Commercial licensing** is available through cuinnovations@cuanschutz.edu (mention ID: CU6153H). Happy to help with integration into your pipeline.

For more information, feel free to `contact me<https://www.jenniferfouquier.com/contact>`_.

Install
========

Source code

Run the command:

`git clone git@github.com:JTFouquier/explana.git`

Next, navigate to directory with `Snakefile` and work from there.



Create Conda environment using Mamba
========================================

You must be inside the directory with `Snakefile` to create the `explana` environment and run your analyses.

You will need to install Mambaforge `download Mambaforge here<https://github.com/conda-forge/miniforge#mambaforge>`_

`mamba env create -f conda_envs/environment.yaml -n explana`

Activate the `explana` conda environment. This is needed to ensure specific software versions are used.

`conda activate explana`


Test install
============

Run EXPLANA with small install demo. If possible, use more cores (I use 6).

`snakemake --configfile config/config-install-test.yaml --cores 2 --use-conda`

Open html report (double click from within directory)

`open workflow-results/EXPLANA-install-test/EXPLANA-report.html`

To retest, remove output folder.

To remove and retest (`rm -r` removes a folder and *cannot* be undone):

`rm -r workflow-results/EXPLANA-install-test/`

If workflow partially runs, rerun `snakemake` command and Snakemake will resume analysis. For major changes, you should rerun the whole workflow.