import os
from Bio import SeqIO
import collections
import sys
import argparse

parser = argparse.ArgumentParser(description='Gathering and filtering files.')
parser.add_argument("--blast_file", action="store", type=str, dest="blast")
parser.add_argument("--reference_file", action="store", type=str, dest="reference_file")
parser.add_argument("--reads", action="store", type=str, dest="reads")
# parser.add_argument("--outreads", action="store", type=str, dest="outreads")
parser.add_argument("--output_dir", action="store", type=str, dest="output_dir")
parser.add_argument("--summary", action="store", type=str, dest="summary")
parser.add_argument("--sample", action="store", type=str, dest="sample")
parser.add_argument("--primers", action="store", type=str, dest="primers")

args = parser.parse_args()

blast_file = str(args.blast)
refs=str(args.reference_file)
reads = str(args.reads)

outdir = str(args.output_dir)
summary = str(args.summary)
barcode=str(args.sample)
primers= str(args.primers)

print(barcode)
hits = collections.defaultdict(list)
with open(blast_file, "r") as f:
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
mod_counter = collections.Counter()
for record in SeqIO.parse(reads,"fastq"):
    try:
        gene,mod= top_gene[record.id].split('_')
        mod_counter[top_gene[record.id]]+=1
    except:
        gene="none"
        mod_counter["none"]+=1
    records[gene].append(record)

seq_dict = {}
for record in SeqIO.parse(refs,"fasta"):
    print(record.id)
    seq_dict[record.id]=record.seq.upper()

primer_coords={}
with open(primers,"r") as f:
    for l in f:
        l=l.rstrip('\n')
        tokens=l.split(',')
        primer_coords[tokens[0]]= (len(tokens[1]), len(tokens[2]))

with open(summary,"w") as fwsum:
    for gene in records:
        if gene != "none":
            fwref = open("pipeline_output/binned/"+barcode+"_bin/primer-schemes/minion/V_"+gene+"/minion.reference.fasta","w")
            fwbed = open("pipeline_output/binned/"+barcode+"_bin/primer-schemes/minion/V_"+gene+"/minion.scheme.bed","w")
            unmodified_count = mod_counter[gene+"_unmodified"]
            modified_count = mod_counter[gene+"_modified"]
            # if unmodified_count > modified_count:
            #     #use the unmodified gene as reference for mapping
            #     sequence= seq_dict[gene+"_unmodified"]
            #     fwref.write(">{}\n{}\n".format(gene,sequence))
            #     fwbed.write("{}\t{}\t{}\t{}_LEFT\t1\n".format(gene,0,30,gene))
            #     fwbed.write("{}\t{}\t{}\t{}_RIGHT\t1\n".format(gene,len(sequence)-30,len(sequence),gene))
            # else:
            sequence= seq_dict[gene+"_modified"]
            sequence= str(sequence).replace("CG","TG")
            fwref.write(">{}\n{}\n".format(gene,sequence))
            fwbed.write("{}\t{}\t{}\t{}_LEFT\t1\n".format(gene,0,primer_coords[gene][0],gene))
            fwbed.write("{}\t{}\t{}\t{}_RIGHT\t1\n".format(gene,len(sequence)-primer_coords[gene][1],len(sequence),gene))

            print("In file:{}\n For gene: {}\n\tUnmodified hit count: {}\n\tModified hit count: {}\n".format(blast_file, gene, unmodified_count, modified_count))
            out_file = outdir + '/' + gene + '.fastq'
            fwsum.write("{}: {} records written to {}. \n".format(gene,len(records[gene]),out_file))
            with open(out_file, "w") as fwseq:
                SeqIO.write(records[gene], fwseq, "fastq")
        fwref.close()
        fwbed.close()

