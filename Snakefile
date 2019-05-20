import os
import collections
from Bio import SeqIO
import pysam

configfile: "config.yaml"
run_name = str(config["run_name"])


##### Target rules #####

rule all:
    input:
        expand("pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.consensus.fasta",gene=config["genes"],barcode=config["barcodes"]),
        expand("pipeline_output/binned/{barcode}_bin/reads/{gene}.fastq", gene=config["genes"],barcode=config["barcodes"])

##### Modules #####
include: "rules/gather.smk"
include: "rules/nanopolish_index.smk"
include: "rules/demultiplex.smk"
include: "rules/bin.smk"
include: "rules/minion.smk"


onstart:
    print("Setting up the artic package")
    shell("cd fieldbioinformatics && python setup.py install")
    shell("export PATH=$PATH:`pwd`/artic")
    shell("cd .. ")

