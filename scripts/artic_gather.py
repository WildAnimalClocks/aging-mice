#adapted from Nick Loman's script gather.py
import sys
from Bio import SeqIO
import tempfile
import os
import glob
import shutil
import pandas as pd
from collections import defaultdict

directory = str(snakemake.params.path_to_fastq)
run_name=str(snakemake.params.run_name)
min_length=int(snakemake.params.min_length)
max_length=int(snakemake.params.max_length)

all_fastq_outfn = "{}_all.fastq".format(run_name)
print(all_fastq_outfn)
all_fastq_outfh = open(all_fastq_outfn, "w")

fastq = defaultdict(list)
for root, dirs, files in os.walk(directory):
    paths = os.path.split(root)
    local_dir = paths[-1]
    fastq[local_dir].extend([root+'/'+f for f in files if f.endswith('.fastq')])

for local_dir, fastq in list(fastq.items()):
    if len(fastq):
        fastq_outfn = "{}_{}.fastq".format(run_name, local_dir)
        outfh = open(fastq_outfn, "w")
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
                    SeqIO.write([rec], outfh, "fastq")
                    SeqIO.write([rec], all_fastq_outfh, "fastq")

                    dups.add(rec.id)
                    uniq += 1

        outfh.close()

        print("%s\t%s\t%s" % (fastq_outfn, total, uniq))

all_fastq_outfh.close()

print("Collecting summary files\n", file=sys.stderr)

dfs = []

summary_outfn = "%s_sequencing_summary.txt" % (run_name)
summaryfh = open(summary_outfn, "w")

for r, d, f in os.walk(directory):
    for fn in f:
        if fn.startswith("sequencing_summary") and fn.endswith(".txt"):
            
            summ_file = r + '/' + fn
            df = pd.read_csv(summaryfn, sep="\t")
            dfs.append(df)

pd.concat(dfs).to_csv(summaryfh, sep="\t", index=False)
summaryfh.close()


