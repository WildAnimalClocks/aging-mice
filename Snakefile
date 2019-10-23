import os
import collections
from Bio import SeqIO
import pysam

configfile: "config.yaml"
run_name = str(config["run_name"])


##### Configuration #####

# trim trailing slashes from paths to avoid snakemake complaining of double '/' in paths
config["output_path"] = config["output_path"].rstrip("/")
config["input_path"] = config["input_path"].rstrip("/")

# todo - check that 'barcode_set' is one of 'native', 'rapid', `pcr` or 'all' and throw error if not
barcode_set = " --native_barcodes"
if str(config["barcode_set"]).lower()=="native":
    barcode_set= " --native_barcodes"
elif str(config["barcode_set"]).lower()=="rapid":
    barcode_set= " --rapid_barcodes"
elif str(config["barcode_set"]).lower()=="pcr":
    barcode_set=" --pcr_barcodes"
elif str(config["barcode_set"]).lower()=="all":
    barcode_set=""

# todo - check that the barcode labels fit the required format (text followed by an integer)
limit_barcodes_to = ""
try:
    barcodes = config["limit_barcodes_to"].split(',')
    barcode_string = ""
    for i in barcodes:
        barcode_string+=" {}".format(i.lstrip("NB").lstrip("BC").lstrip("barcode"))
    limit_barcodes_to = " --limit_barcodes_to {}".format(barcode_string)
except:
    pass
    # print("No barcodes given, will search for all barcodes within the {} barcode list.".format(config["barcode_set"]))


discard_unassigned = ""
if str(config["discard_unassigned"]).lower()=="true":
    discard_unassigned= " --discard_unassigned"
else:
    discard_unassigned= ""

# todo - check that the value is true or false
require_two_barcodes = ""
if str(config["require_two_barcodes"]).lower()=="false":
    require_two_barcodes= ""
else:
    require_two_barcodes= " --require_two_barcodes"

# todo - check that the value is true or false
split_reads = ""
if str(config["split_reads"]).lower()=="true":
    split_reads= " --no_split"
else:
    split_reads= " --no_split"

# todo - check that the value is true or false
discard_middle = ""
if str(config["discard_middle"]).lower()=="true":
    discard_middle= " --discard_middle"
else:
    discard_middle= ""



##### Target rules #####

rule all:
    input:
        "pipeline_output/cpg_report.csv",
        config["output_path"]+ "/{}.csv".format(run_name),
        expand(config["output_path"]+ "/{barcode}_bin/{barcode}.fastq", barcode=config["barcodes"]),
        expand(config["output_path"]+ "/{barcode}_bin/reads/{gene}.fastq", gene=config["genes"], barcode=config["barcodes"])


##### Modules #####
include: "rules/gather.smk"
include: "rules/nanopolish_index.smk"
include: "rules/demultiplex.smk"
include: "rules/bin.smk"
include: "rules/minion.smk"
include: "rules/count.smk"

onstart:
    print("Setting up the artic package")
    shell("cd fieldbioinformatics && python setup.py install")
    shell("export PATH=$PATH:`pwd`/artic")
    shell("cd .. ")

