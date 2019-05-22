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

<img src="https://github.com/aineniamh/aging-mice/blob/master/one_sample_dag.svg">


1. setup ``artic fieldbioinformatics`` package \
Automatic setup of this on startup of the pipeline. Gives the user access to ``artic minion`` script for step below.
2. gather \
Parses all of the basecalled fastq files from ``guppy``, applies a length filter that can be customised in the ``config.yaml`` file and writes the reads to a single file ``run_name_all.fastq``. This script also searches the fastq directories for ``sequencing_summary`` files and combines them into a single file: ``run_name_sequencing_summary.txt``. These files will be output in the ``pipeline_output`` directory.
3. demultiplex_qcat \
For each read in the ``run_name_all.fastq`` file, identifies barcodes and outputs reads into respective files, binned by barcode. These files appear in the ``demultiplexed`` directory, in ``pipeline_output``.
4. blastn \
For each ``barcode.fastq`` file, each read is blasted against a database containing the 5 genes of interest in modified and unmodified form.
5. bin \
This step parses each blast output and assesses for each read what the best blast hit is. The reads are then binned by gene and a count of 'modified vs unmodified' best blast hits for each barcode for each gene is performed. This determines which reference (modified or unmodified) is most suited to take forward into nanopolish for each barcode for each gene. It also creates the respective bed file for the ``artic minion`` pipeline to use.
6. nanopolish_index \
Creates the nanopolish index necessary for running nanopolish in the next step. It accesses the gathered fastq and sequencing summary files from step 2 and also the signal-level fast5 data.
7. artic_minion \
The ``artic minion`` pipeline, written by Nick Loman, is then run for each ``barcode_gene`` combination in order to generate a high-quality consensus sequence, using an approach informed by signal-level data. This pipeline performs the following steps:
    * Maps against a given reference and sorts reads using ``bwa`` and ``samtools`` respectively.
    * Runs the ``artic align_trim`` script. This script takes in a bed file and your alignment and assesses whether the primers are correctly paired according to the bed file, discarding reads that are not, and normalises the read coverage across the genome. It is run twice, first to trim off the barcodes and the primers and second to just trim off the barcodes.
    * Loads the ``nanopolish index`` created in step 8.
    * Runs ``nanopolish variants`` twice, on the barcode-and-primer-trimmed bam and on the barcode-trimmed bam.
    * Generates a variant frequency plot.
    * Runs ``margin_cons``, a custom script that filters the variants, masking sites that do not reach the depth threshold of 20 and do not reach a quality threshold of 200, and produces a consensus sequence with 'N' masking on the relevant sites. It uses the vcf from nanopolish without primer-trimming but the primer-trimmed bam file so that primer sequences do not count towards depth calculation. A report is also generated.
8. organise_minion_output \
Moves the output for each gene into the relevant barcode directory ``pipeline_output/minion_output/barcodeX_bin/``.
9. count \
Parses the output vcf files and produces a report, ``cpg_report.csv`` in the ``pipeline_output`` directory, with information about percentage reads and read counts of modified/ unmodified CpG sites.
