# aging-mice

Repository for a snakemake pipeline for running analysis of aging-mice MinION sequencing project.


## setup

Although not a requirement, an install of conda will make the setup of this pipeline on your local machine much more streamlined. I have created an ``aging-mice`` conda environment which will allow you to access all the software required for the pipeline to run. To install conda, visit here https://conda.io/docs/user-guide/install/ in a browser. 

> *Recommendation:* Install the `64-bit Python 3.6` version of Miniconda

Once you have a version of conda installed on your machine, clone this repository by typing into the command line:

```bash
git clone https://github.com/aineniamh/aging-mice.git
```

Build the conda environment by typing:

```bash
conda env create -f aging-mice/environment.yml
```

To activate the environment, type:

```bash
source activate aging-mice
```

To deactivate the environment, enter:

```bash
conda deactivate
```

## customising the pipeline

Inside the ``config.yaml`` file, you can insert your MinION ``run_name`` and the barcode and gene names. I have filled in the 5 gene names from the original PCR and MinION run, but this can be modified here in the future if the study design changes.