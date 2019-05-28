import os
import numpy as np
from Bio import AlignIO
import collections

fw = open("references/cpg_sites.csv","w")
fw.write("gene,mod_id,position(1-based)\n")
mod_c = 0
gene_mod = 0
for record in SeqIO.parse("references/genes.all.fasta","fasta"):
    if record.id.endswith("unmodified"):
        print(record.id, record.seq.upper())
        gene = record.id.rstrip("_unmodified")
        seq = record.seq.upper()
        print(primer_coords[gene])
        for i in range(len(seq)):
#             print(i)
            if i > primer_coords[gene][0] and i < len(seq)-primer_coords[gene][1]:
                if seq[i:i+2]=="CG":
                    mod_c +=1
                    print(i, seq[i:i+2])
                    fw.write("{},mod_{},{}\n".format(gene,mod_c,i+1))
fw.close()
