# aging-mice

Repository for a snakemake pipeline for running analysis for the aging-mice MinION sequencing project.

## background

 DNA was extracted from 24 mice, of varying ages and bisulfite treated. PCR was run in order to amplify a number of targeted genes using primers that target modified sequence. Sequencing was performed using the MinION and the fast5 reads were basecalled with MinKNOW. This pipeline starts with the basecalled reads, gathers them together and applies a read length filter. It then demultiplexes them using ``porechop`` with strict double barcoding. ``ParaMethER`` (Parasail Methylation Estimator for Regression) then uses parasail to align each read against a panel of references and it selects the best reference from the set for each read (applying a minimal coverage (80%) and identity (55%) threshold). ``ParaMethER`` then extracts the known CpG positions from the alignment and counts them up for each sample. It also assesses a background error rate of T->C misassignment and includes this in the report, per gene, per sample. Reports are output into the output directory specified in the config file.

## setup

Although not a requirement, an install of conda will make the setup of this pipeline on your local machine much more streamlined. I have created an ``epi-clock`` conda environment which will allow you to access all the software required for the pipeline to run. To install conda, visit here https://conda.io/docs/user-guide/install/ in a browser. 

> *Recommendation:* Install the `64-bit Python 3.6` version of Miniconda

Once you have a version of conda installed on your machine, clone this repository by typing into the command line:

```bash
git clone https://github.com/aineniamh/aging-mice.git
```

Build the conda environment by typing:

```bash
conda env create -f environment.yaml
```

To activate the environment, type:

```bash
source activate epi-clock
```

To deactivate the environment (after you have finished your analysis), enter:

```bash
conda deactivate
```

## customising the pipeline

To run the analysis using snakemake, you will need to customise the ``config.yaml`` file.

Inside the ``config.yaml`` file, you can change your MinION ``run_name`` and the barcode and gene names. I have filled in the 4 gene names from the original PCR and MinION run, but this can be modified here in the future if the study design changes.

Ensure ``input_path`` points to where your fastq reads are.

## running the pipeline

To start the pipeline, in a terminal window in the aging-mice directory, simply enter:

```bash
snakemake
```

If you wish to run your pipeline using more than one core (**recommended**), enter:

```bash
snakemake --cores X
```

where X is the number of threads you wish to run.
