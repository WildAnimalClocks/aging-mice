import os
from Bio import SeqIO
import collections
import sys
import argparse

parser = argparse.ArgumentParser(description='Gathering and filtering files.')
parser.add_argument("--blast_file", action="store", type=str, dest="blast")
parser.add_argument("--reference_file", action="store", type=str, dest="reference_file")
parser.add_argument("--reads", action="store", type=str, dest="reads")
parser.add_argument("--output_dir", action="store", type=str, dest="output_dir")
parser.add_argument("--summary", action="store", type=str, dest="summary")
parser.add_argument("--sample", action="store", type=str, dest="sample")

args = parser.parse_args()

outdir = str(args.output_dir)
summary = str(args.summary)
barcode=str(args.sample)

print(barcode)
hits = collections.defaultdict(list)
with open(str(args.blast_file), "r") as f:
    for l in f:
        tokens= l.rstrip('\n').split(',')
        read = tokens[0]
        hit = tokens[1]
        score= tokens[-1] 
        hits[read].append((hit,score))

top_gene = {}
for i in hits:
    top_hit = sorted(hits[i], key = lambda x : float(x[1]), reverse = True)[0]
    top_gene[i]=top_hit[0]
hits.clear()

records = collections.defaultdict(list)
for record in SeqIO.parse(str(args.reads),"fastq"):
    try:
        gene= top_gene[record.id].split('_')[0]
    except:
        gene="none"
    records[gene].append(record)

seq_dict = {}
for record in SeqIO.parse(args.reference_file,"fasta"):
    print(record.id)
    seq_dict[record.id]=record.seq.upper()

primer_coords={}
with open(primers,"r") as f:
    for l in f:
        l=l.rstrip('\n')
        tokens=l.split(',')
        primer_coords[tokens[0]]= (len(tokens[1]), len(tokens[2]))

with open(summary,"w") as fwsum:
    fwsum.write("Gene,Num_reads\n")
    for gene in records:
        fwsum.write("{},{}\n".format(gene,len(records[gene])))
        print("In file:{}\n For gene: {}\n\tNumber of reads: {}\n".format(args.blast_file, gene, len(records[gene])))

        if gene != "none":

            out_file = outdir + '/' + gene + '.fastq'
            with open(out_file, "w") as fwseq:
                SeqIO.write(records[gene], fwseq, "fastq")
