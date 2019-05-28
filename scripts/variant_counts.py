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
parser.add_argument("--cpg_wide", action="store", type=str, dest="cpg_wide")
parser.add_argument("--ages", action="store", type=str, dest="ages")
args = parser.parse_args()

age_dict = {}
with open(str(args.ages),"r") as f:
    for l in f:
        l = l.rstrip('\n')
        tokens = l.split(',')
        age_dict[tokens[0]]=tokens[1]

mod_list= []
mod_string = 'barcode,age,'
with open(str(args.cpg),"r") as f:
    for l in f:
        if not l.startswith("gene"):
            l = l.rstrip('\n')
            tokens = l.split(',')
            mod_list.append(tokens[0]+"_"+tokens[2])
            mod_string+=tokens[0]+"_"+tokens[2]+","
mod_string=mod_string.rstrip(',')
mod_string += "\n"
    
file_mods=defaultdict(dict)
for r,d,f in os.walk("pipeline_output_2/minion_output"):
    for filename in f:
        if not filename.endswith("primertrimmed.vcf") and filename.endswith(".vcf"):
            print("Counting CpG modifications in file: {}".format(filename))

            barcode,gene=filename.rstrip('.vcf').split('_')
            print(barcode,gene)
            with open(r+ '/' + filename,"r") as f:
                for l in f:
                    l=l.rstrip('\n')
                    if l.startswith('#'):
                        pass
                    else:
                        tokens=l.split()
                        chrom,pos,dot,ref,alt,qual,pass_info,info,fmat,no=tokens
                        if ref=="T" and alt=="C":
                            info_list=info.split(';')
                            num_var_reads=int(info_list[0].lstrip("BaseCalledReadsWithVariant="))
                            supp_fraction=float(info_list[4].lstrip("SupportFraction="))
                            total_reads=int(info_list[2].lstrip("TotalReads="))
                            print(gene, pos)
                            file_mods[(barcode,age_dict[barcode])][gene+'_'+pos]=(supp_fraction,num_var_reads,total_reads)

fwide = open(str(args.cpg_wide),"w")
fwide.write(mod_string)

fw = open(str(args.out_file),"w")
fw.write("barcode,age,gene_position,support_fraction,num_reads_c,total_reads\n")

for key in file_mods:

    barcode,age=key
    fwide.write(barcode+','+age+',')
    l = ''
    for i in mod_list:
        try:
            supp_fraction=str(file_mods[key][i][0])
            l+=supp_fraction+','
            fw.write("{},{},{},{},{},{}\n".format(barcode,age,i,file_mods[key][i][0],file_mods[key][i][0],file_mods[key][i][0]))
        except:
            l+='1,'
            fw.write("{},{},{},{},{},{}\n".format(barcode,age,i,1,"NA","NA"))
            
    l=l.rstrip(',')
    fwide.write(l+'\n')
fwide.close()
fw.close()