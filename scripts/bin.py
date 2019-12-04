import os
from Bio import SeqIO
import collections
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Parsing blast output.')

    parser.add_argument("--blast_file", action="store", type=str, dest="blast")
    parser.add_argument("--references", action="store", type=str, dest="references")

    parser.add_argument("--reads", action="store", type=str, dest="reads")

    parser.add_argument("--output_dir", action="store", type=str, dest="output_dir")
    parser.add_argument("--report", action="store", type=str, dest="report")

    parser.add_argument("--barcode", action="store", type=str, dest="barcode")

    return parser.parse_args()

def load_references(refs):
    ref_dict = {}
    for record in SeqIO.parse(refs, "fasta"):
        ref_dict[record.id]= []

    return ref_dict

def get_hits(blast):
    hits = collections.defaultdict(list)
    with open(str(blast), "r") as f:
        for l in f:
            tokens= l.rstrip('\n').split(',')
            read = tokens[0]
            hit = tokens[1]
            score= tokens[-1] 
            hits[read].append((hit,score))
    return hits

def greatest_hits(blast, ref_dict):
    counts = {
        "total": 0
    }
    hits = get_hits(blast)
    for i in hits:
        counts["total"] +=1
        top_hit = sorted(hits[i], key = lambda x : float(x[1]), reverse = True)[0]
        ref_dict[top_hit[0]].append(i)
    return counts
    
def get_records(read_names, reads):
    records = []
    for record in SeqIO.parse(reads, "fastq"):
        if record.id in read_names:
            records.append(record)

    return records

def get_unmapped(total_hits, reads):
    total = 0
    for record in SeqIO.parse(reads, "fastq"):
        total +=1
    return total-total_hits

if __name__ == '__main__':

    args = parse_args()

    ref_dict = load_references(args.references)

    counts = greatest_hits(args.blast, ref_dict)
    unmapped = get_unmapped(counts["total"], str(args.reads))

    freport = open(str(args.report), "w")
    freport.write(f"{args.barcode},unmapped,{unmapped}\n")

    for ref in ref_dict:
        freport.write(f"{args.barcode},{ref},{len(ref_dict[ref])}\n")
        with open(str(args.output_dir)+'/'+ref+".fastq", "w") as fw:
            
            gene_records = get_records(ref_dict[ref], str(args.reads))
            SeqIO.write(gene_records, fw, "fastq")

    freport.close()


# primer_coords={}
# with open(str(args.primers),"r") as f:
#     for l in f:
#         l=l.rstrip('\n')
#         tokens=l.split(',')
#         primer_coords[tokens[0]]= (len(tokens[1]), len(tokens[2]))

# with open(report,"w") as fwsum:
#     fwsum.write("Gene,Num_reads\n")
#     for gene in records:
#         fwsum.write("{},{}\n".format(gene,len(records[gene])))
#         print("In file:{}\n For gene: {}\n\tNumber of reads: {}\n".format(args.blast_file, gene, len(records[gene])))

#         if gene != "none":
#             with open(outdir + '/' + gene + ".primers.bed","w") as fw:
#                 fw.write("{},0,{}\n".format(gene, primer_coords[gene][0]))
#                 fw.write("{},{},{}\n".format(gene, (len(seq_dict[gene])-primer_coords[gene][1]), len(seq_dict[gene])))

#             with open(outdir + '/' + gene + ".fasta","w") as fw:
#                 fw.write(">{}\n{}\n".format(gene, seq_dict[gene]))

#             with open(outdir + '/' + gene + '.fastq', "w") as fwseq:
#                 SeqIO.write(records[gene], fwseq, "fastq")
