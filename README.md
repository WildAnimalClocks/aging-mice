# aging-mice

Repository for a snakemake pipeline for running analysis for the aging-mice MinION sequencing project.

## background

 DNA was extracted from 24 mice, of varying ages and bisulfite treated. PCR was run in order to amplify a number of targeted genes using primers that target modified sequence. Sequencing was performed using the MinION and the fast5 reads were basecalled with MinKNOW. This pipeline starts with the basecalled reads, gathers them together and applies a read length filter. It then demultiplexes them using ``qcat`` and bins them by a defined list of genes using ``BLAST`` and a custom python script. From there, the reads are mapped, sorted, indexed and variant called using a combination of samtools, custom scripts and nanopolish. Custom scripts then collect the information from the variant calls and output a summary csv.

## setup

Although not a requirement, an install of conda will make the setup of this pipeline on your local machine much more streamlined. I have created an ``aging-mice`` conda environment which will allow you to access all the software required for the pipeline to run. To install conda, visit here https://conda.io/docs/user-guide/install/ in a browser. 

> *Recommendation:* Install the `64-bit Python 3.6` version of Miniconda

Once you have a version of conda installed on your machine, clone this repository by typing into the command line:

```bash
git clone https://github.com/aineniamh/aging-mice.git
```

Build the conda environment by typing:

```bash
conda env create -f aging-mice/envs/aging-mice.yaml
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

To run the analysis using snakemake, you may need to customise the ``config.yaml`` file.

Inside the ``config.yaml`` file, you can change your MinION ``run_name`` and the barcode and gene names. I have filled in the 5 gene names from the original PCR and MinION run, but this can be modified here in the future if the study design changes.

Ensure ```path_to_fast5``` and ```path_to_fastq``` point to where your data is.

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

# pipeline description

1. gather \
Parses all of the basecalled fastq files from ``guppy``, applies a length filter that can be customised in the ``config.yaml`` file and writes the reads to a single file ``run_name_all.fastq``. This script also searches the fastq directories for ``sequencing_summary`` files and combines them into a single file: ``run_name_sequencing_summary.txt``. These files will be output in the ``pipeline_output`` directory.
2. demultiplex_qcat \
For each read in the ``run_name_all.fastq`` file, identifies barcodes and outputs reads into respective files, binned by barcode. These files appear in the ``demultiplexed`` directory, in ``pipeline_output``.
3. blastn \
For each ``barcode.fastq`` file, each read is blasted against a database containing the 5 genes of interest in modified and unmodified form.
4. bin \
This step parses each blast output and assesses for each read what the best blast hit is. The reads are then binned by gene and a count of 'modified vs unmodified' best blast hits for each barcode for each gene is performed. This determines which reference (modified or unmodified) is most suited to take forward into nanopolish for each barcode for each gene.
5. minimap2, sort the reads and index them using ``samtools``.
6. Downsample and filter by read quality using a custom script (still need to write this).
7. ``nanopolish index`` and run ``nanopolish``.