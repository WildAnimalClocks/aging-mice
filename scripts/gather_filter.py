#adapted from Nick Loman's script gather.py
import sys
from Bio import SeqIO
import tempfile
import os
import pandas as pd
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='Gathering and filtering files.')
parser.add_argument("--min_length", action="store", type=int, dest="min_length")
parser.add_argument("--max_length", action="store", type=int, dest="max_length")
parser.add_argument("--run_name", action="store", type=str, dest="run_name")
parser.add_argument("--path_to_fastq", action="store", type=str, dest="path_to_fastq")
parser.add_argument("--output_directory", action="store", type=str, dest="output_directory")

params = parser.parse_args()

directory = params.path_to_fastq

min_length=params.min_length
max_length=params.max_length


all_fastq_outfn = "{}/{}.fastq".format(params.output_directory, params.run_name)
all_fastq_outfh = open(all_fastq_outfn, "w")

fastq = defaultdict(list)
for root, dirs, files in os.walk(directory):
    paths = os.path.split(root)
    local_dir = paths[-1]
    fastq[local_dir].extend([root+'/'+f for f in files if f.endswith('.fastq')])

for local_dir, fastq in list(fastq.items()):
    if len(fastq):
        print("Processing {} files in {}".format(len(fastq), local_dir))

        dups = set()
        uniq = 0
        total = 0    
        for f in fastq:
            for rec in SeqIO.parse(open(f), "fastq"):
                if len(rec) > max_length:
                    continue
                if len(rec) < min_length:
                    continue

                total += 1
                if rec.id not in dups:
                    SeqIO.write([rec], all_fastq_outfh, "fastq")

                    dups.add(rec.id)
                    uniq += 1

        print("Processed: {}\t{}".format(total, uniq))

all_fastq_outfh.close()



