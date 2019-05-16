import os
import collections
from Bio import SeqIO

configfile: "config.yaml"
run_name = str(config["run_name"])


##### Target rules #####

rule all:
    input:
        expand("pipeline_output/binned/{barcode}/binning_report.txt",barcode=config["barcodes"]),
        expand("pipeline_output/binned/{barcode}/mapped/{gene}.sam",gene=config["genes"],barcode=config["barcodes"])

##### Modules #####
include: "rules/gather.smk"
include: "rules/demultiplex.smk"
include: "rules/mapping.smk"
include: "rules/nanopolish_index.smk"
include: "rules/bin.smk"
