import os


configfile: "config.yaml"
run_name = str(config["run_name"])


##### Target rules #####

rule all:
    input:
        expand("demultiplexed/{barcode}.fastq",barcode=config["barcodes"])


##### Modules #####
include: "rules/gather.smk"
include: "rules/demultiplex.smk"
#include: "rules/call.smk"
#include: "rules/bin.smk"
