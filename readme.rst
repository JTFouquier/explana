**© 2025 Jennifer Fouquier LLC. ALL RIGHTS RESERVED**

EXPLANA: A user-friendly workflow for EXPLoratory ANAlysis and feature selection in cross-sectional and longitudinal studies
=====================================================================================================================================

.. image:: https://github.com/JTFouquier/explana/blob/main/images/report-screenshot.png
   :width: 50%

EXPLANA was developed to streamline identification of features that relate to outcome variables in longitudinal, and cross-sectional, microbiome studies. It has broad application to different types of data.

Outcome variables can be numerical or categorical (binary - yes/no) and input features can be numerical or categorical. If order of categorical value changes impact a response, this can be identified.

Features are ranked by their importance, and feature impact on response is provided. An interactive .html report is generated to visually and textually summarize results.

This is an exploratory method and this should be made clear when communicating results.

User Guide
===========

`explana.io/documentation <https://www.explana.io/documentation>`_

Licensing/Use
=============

Software can be installed for academic/non-profit use.

For commercial licensing or consulting, email jennifer@explana.io or `schedule time to chat <https://www.jenniferfouquier.com/booking-calendar/availability>`_.

Install
========

Source code

Run the command:

``git clone git@github.com:JTFouquier/explana.git``

Next, navigate to directory with ``Snakefile`` and work from there.

Create Conda environment using Mamba
========================================

You must be inside the directory with ``Snakefile`` to create the ``explana`` environment and run your analyses.

You will need to install Mambaforge `download Mambaforge here <https://github.com/conda-forge/miniforge#mambaforge>`_

``mamba env create -f conda_envs/environment.yaml -n explana``

Activate the ``explana`` conda environment. This is needed to ensure specific software versions are used.

``conda activate explana``

Test install
============

Run EXPLANA with small install demo. If possible, use more cores.

``snakemake --configfile config/config-install-test.yaml --cores 2 --use-conda``

Open html report (double click from within directory)

``open workflow-results/EXPLANA-install-test/EXPLANA-report.html``

To retest, remove output folder.

To remove and retest (``rm -r`` removes a folder and *cannot* be undone):

``rm -r workflow-results/EXPLANA-install-test/``

If workflow partially runs, rerun ``snakemake`` command and Snakemake will resume analysis. For major changes, you should rerun the whole workflow.