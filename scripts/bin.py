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
parser.add_argument("--output-dir", action="store", type=str, dest="output_dir")
parser.add_argument("--summary", action="store", type=str, dest="summary")

args = parser.parse_args()
if not os.path.exists("pipeline_output/binned"):
    os.mkdir("pipeline_output/binned")

if not os.path.exists("pipeline_output/binned/{barcode}_bin/"):
    os.mkdir("pipeline_output/binned/{barcode}_bin/")
    os.mkdir("pipeline_output/binned/{barcode}_bin/primer-schemes")
    os.mkdir("pipeline_output/binned/{barcode}_bin/primer-schemes/minion")

blast_file = str(args.blast)
refs=str(args.reference_file)
reads = str(args.reads)
outdir = str(args.output_dir)
summary = str(args.summary)

hits = collections.defaultdict(list)
with open(blast_file, "r") as f:
    for l in f:
        tokens= l.rstrip('\n').split(',')
        read = tokens[0]
        hit = tokens[1]
        score= tokens[-1] #check these
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
    seq_dict[record.id]=record.seq.upper()

with open(summary,"w") as f:
    for gene in records:
        if gene != "none":  
            os.mkdir("pipeline_output/binned/{barcode}_bin/primer-schemes/minion/V_"+gene)

            unmodified_count = mod_counter[gene+"_unmodified"]
            modified_count = mod_counter[gene+"_modified"]
            if unmodified_count >= modified_count:
                #use the unmodified gene as reference for mapping
                sequence= seq_dict[gene+"_unmodified"]
                with open("pipeline_output/binned/{barcode}_bin/primer-schemes/minion/V_"+gene+"/minion.reference.fasta","w") as fnano:
                    fnano.write(">{}_unmodified\n{}\n".format(gene,sequence))
                with open("pipeline_output/binned/{barcode}_bin/primer-schemes/minion/V_"+gene+"/minion.scheme.bed","w") as fnano:
                    fnano.write("{}\t{}\t{}\t{}_LEFT\t1\n".format(gene,0,30,gene))
                    fnano.write("{}\t{}\t{}\t{}_RIGHT\t1\n".format(gene,len(sequence)-30,len(sequence),gene))
            else:
                sequence= seq_dict[gene+"_modified"]
                with open("pipeline_output/binned/{barcode}_bin/primer-schemes/minion/V_"+gene+"/minion.reference.fasta","w") as fnano:
                    fnano.write(">{}_modified\n{}\n".format(gene,sequence))
                with open("pipeline_output/binned/{barcode}_bin/primer-schemes/minion/V_"+gene+"/minion.scheme.bed","w") as fnano:
                    fnano.write("{}\t{}\t{}\t{}_LEFT\t1\n".format(gene,0,30,gene))
                    fnano.write("{}\t{}\t{}\t{}_RIGHT\t1\n".format(gene,len(sequence)-30,len(sequence),gene))

            print("In file:{}\n For gene: {}\n\tUnmodified hit count: {}\n\tModified hit count: {}\n".format(blast_file, gene, unmodified_count, modified_count))
            out_file = outdir + '/' + gene + '.fastq'
            f.write("{}: {} records written to {}. \n".format(gene,len(records[gene]),out_file))
            with open(out_file, "w") as fw:
                SeqIO.write(records[gene], fw, "fastq")