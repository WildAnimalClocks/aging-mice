import sys
from Bio import SeqIO
import tempfile
import os
import pandas as pd
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='Counting modifications.')
parser.add_argument("--out_file", action="store", type=str, dest="out_file")
parser.add_argument("--cpg_info", action="store", type=str, dest="cpg")
args = parser.parse_args()

mod_dict=defaultdict(list)
with open(str(args.cpg),"r") as f:
    for l in f:
        if not l.startswith("gene"):
            l = l.rstrip('\n')
            tokens = l.split(',')
            mod_dict[tokens[0]].append(tokens[3])
    
fw = open(str(args.out_file),"w")
fw.write("barcode,gene,position,pcent_c,num_reads_c,total_reads\n")

for r,d,f in os.walk("pipeline_output/minion_output"):
    for filename in f:
        if not filename.endswith("primertrimmed.vcf") and filename.endswith(".vcf"):
            print("Counting CpG modifications in file: {}".format(filename))

            barcode,gene=filename.rstrip('.vcf').split('_')
            
            file_mods = {}
            for i in mod_dict[gene]:
                if int(i) > 30:
                    file_mods[i]=(100,200,200)
            
            with open(r+ '/' + filename,"r") as f:
                for l in f:
                    l=l.rstrip('\n')
                    if l.startswith('#'):
                        pass
                    else:
                        tokens=l.split()
                        chrom,pos,dot,ref,alt,qual,pass_info,info,fmat,no=tokens
                        info_list=info.split(';')
                        num_var_reads=int(info_list[0].lstrip("BaseCalledReadsWithVariant="))
                        pcent_c=round(100-(100*float(info_list[1].lstrip("BaseCalledFraction="))),4)
                        total_reads=int(info_list[2].lstrip("TotalReads="))

                        file_mods[pos]=(pcent_c,total_reads-num_var_reads,total_reads)
            for i in file_mods:
                mod_info = file_mods[i]  

                fw.write("{},{},{},{},{},{}\n".format(barcode,gene, i, file_mods[i][0],file_mods[i][1],file_mods[i][2]))

fw.close()

