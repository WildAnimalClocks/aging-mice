import os
import numpy as np
from Bio import AlignIO
import collections

fw = open("references/cpg_sites.csv","w")
fw.write("gene,mod_id,mod_along_gene,position\n")
c_all=0
for r,d,f in os.walk("references/mod_and_unmod/"):
    for filename in f:
        print(filename)
        a=AlignIO.read(r+'/'+filename, "fasta")
        print(len(a[0]))
        c=0
        for i in range(len(a[0])):
            base = a[:, i:i+2]
            if base[0].seq.upper()=="CG":
                c+=1
                c_all+=1

                fw.write("{},mod_{},mod_{},{}\n".format(filename.rstrip('.fasta'), c_all,c, i+1))
fw.close()
        