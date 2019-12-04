import argparse
from Bio import SeqIO
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Parse minimap paf file, create gene bins.')

    parser.add_argument("--paf_file", action="store", type=str, dest="paf_file")
    parser.add_argument("--reads", action="store", type=str, dest="reads")

    parser.add_argument("--output_dir", action="store", type=str, dest="output_dir")
    parser.add_argument("--references", action="store", type=str, dest="references")
    parser.add_argument("--report", action="store", type=str, dest="report")
    parser.add_argument("--barcode", action="store", type=str, dest="barcode")

    return parser.parse_args()


def load_references(refs):
    ref_dict = {}
    for record in SeqIO.parse(refs, "fasta"):
        ref_dict[record.id]= []

    return ref_dict

def parse_line(line):

    values = {}
    tokens = line.rstrip('\n').split()
    values["read_name"], values["read_len"] = tokens[:2]    
    values["ref_hit"], values["ref_len"], values["coord_start"], values["coord_end"], values["matches"], values["aln_block_len"] = tokens[5:11]

    return values


def parse_paf(paf, ref_dict):
    counts = {
        "unmapped": 0,
        "split_hit":0,
        "total": 0
    }

    with open(str(paf),"r") as f:
        for line in f:
            counts["total"]+=1
            mapping = parse_line(line)
            if mapping["ref_hit"] == "*":
                counts["unmapped"]+=1
            else:
                if mapping["read_name"] in ref_dict[mapping["ref_hit"]]:
                    counts["split_hit"] +=1
                else:
                    ref_dict[mapping["ref_hit"]].append(mapping["read_name"])

    return counts

def get_records(read_names, reads):
    records = []
    for record in SeqIO.parse(reads, "fastq"):
        if record.id in read_names:
            records.append(record)
    return records

if __name__ == '__main__':

    args = parse_args()

    ref_dict = load_references(args.references)
    
    counts = parse_paf(args.paf_file, ref_dict)
    
    freport = open(str(args.report), "w")
    unmapped = counts["unmapped"]
    freport.write(f"{args.barcode},unmapped,{unmapped}\n")

    for ref in ref_dict:
        freport.write(f"{args.barcode},{ref},{len(ref_dict[ref])}\n")
        with open(str(args.output_dir)+'/'+ref+".fastq", "w") as fw:
            
            gene_records = get_records(ref_dict[ref], str(args.reads))
            SeqIO.write(gene_records, fw, "fastq")
    freport.close()




        
