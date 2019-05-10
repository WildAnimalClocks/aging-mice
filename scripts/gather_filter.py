#adapted from Nick Loman's script gather.py
import sys
from Bio import SeqIO
import tempfile
import os
import pandas as pd
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='Gathering and filtering files.')
parser.add_argument("--min-length", action="store", type=int, dest="min_length")
parser.add_argument("--max-length", action="store", type=int, dest="max_length")
parser.add_argument("--run-name", action="store", type=str, dest="run_name")
parser.add_argument("--path-to-fastq", action="store", type=str, dest="path_to_fastq")

params = parser.parse_args()

directory = params.path_to_fastq
run_name=params.run_name
min_length=params.min_length
max_length=params.max_length

all_fastq_outfn = "{}_all.fastq".format(run_name)
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

print("Collecting summary files\n", file=sys.stderr)

dfs = []

summary_outfn = "{}_sequencing_summary.txt".format(run_name)
summaryfh = open(summary_outfn, "w")

for r, d, f in os.walk(directory):
    for fn in f:
        if fn.endswith(".txt"):
            print(fn)
            summ_file = r + '/' + fn
            df = pd.read_csv(summ_file, sep="\t")
            dfs.append(df)
pd.concat(dfs).to_csv(summaryfh, sep="\t", index=False)
summaryfh.close()


