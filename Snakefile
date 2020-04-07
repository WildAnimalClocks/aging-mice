import os
import collections
from Bio import SeqIO

configfile: "config_shag.yaml"
##### Configuration #####

# trim trailing slashes from paths to avoid snakemake complaining of double '/' in paths
if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip("/")
else:
    config["output_path"] = "analysis"

if config.get("input_path"):
    config["input_path"] = config["input_path"].rstrip("/")


run_name = config["run_name"]

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
        config["output_path"] + "/"+run_name+".fastq",
        expand(config["output_path"]+ "/demultiplexed_reads/{barcode}.fastq", barcode=config["barcodes"]),
        config["output_path"] + "/reports/cpg_counts.csv",
        config["output_path"] +"/reports/cpg_wide.csv"

##### Modules #####
include: "rules/gather.smk"
include: "rules/demultiplex.smk"
include: "rules/count.smk"


