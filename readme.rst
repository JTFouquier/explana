**© 2021-2026 Jennifer Fouquier LLC. ALL RIGHTS RESERVED**

EXPLANA: A user-friendly workflow for EXPLoratory ANAlysis and feature selection in cross-sectional and longitudinal studies
=====================================================================================================================================

.. image:: https://github.com/JTFouquier/explana/blob/main/images/report-screenshot.png
   :width: 50%

EXPLANA was developed to streamline identification of key features that relate to outcomes of interest in cross-sectional and longitudinal microbiome studies. It has broad application to different types of data.

Outcome variables can be numerical or categorical (binary - yes/no) and input features can be numerical or categorical. EXPLANA has novel functionality to detect if categorical order impacts response over time (e.g., green_blue differs from blue_green).

Features are ranked by their importance, and feature impact on response is provided. An interactive .html report is generated to visually and textually summarize results.

This is an exploratory method and this should be made clear when communicating results.

**Peer-Reviewed Publication:**

Jennifer Fouquier, Maggie Stanislawski, John O’Connor, Ashley Scadden, Catherine Lozupone, EXPLANA: a user-friendly workflow for EXPLoratory ANAlysis and feature selection in cross-sectional and longitudinal microbiome studies, Bioinformatics, Volume 42, Issue 1, January 2026, btaf658, `https://doi.org/10.1093/bioinformatics/btaf658 <https://doi.org/10.1093/bioinformatics/btaf658>`_

User Guide
===========

`explana.io/documentation <https://www.explana.io/documentation>`_

Licensing/Use
=============

Software can be installed for academic and non-commercial use. 

Questions about use or installation: jennifer.fouquier@cuanschutz.edu

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
