import os
import collections
from Bio import SeqIO

configfile: "config.yaml"
run_name = str(config["run_name"])


##### Target rules #####

rule all:
    input:
        expand("binned/{barcode}/binning_report.txt",barcode=config["barcodes"])
        #expand("binned/{barcode}/{gene}.fastq",gene=config["genes"],barcode=config["barcodes"])

##### Modules #####
include: "rules/gather.smk"
include: "rules/demultiplex.smk"
#include: "rules/call.smk"
include: "rules/bin.smk"
