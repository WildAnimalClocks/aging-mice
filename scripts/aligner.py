import argparse
import os
import sys

import parasail
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='Parsing blast output.')

    parser.add_argument("--reads", action="store", type=str, dest="reads")
    parser.add_argument("--reference", action="store", type=str, dest="reference")
    parser.add_argument("--cpg_sites", action="store", type=str, dest="cpg")
    

    return parser.parse_args()

def parse_cigar(cigar, reference):

    """
        #M   alignment match (can be a sequence match or mismatch)
        #I   insertion to the reference
        #D   deletion from the reference
        #N   skipped region from the reference
        #S   soft clipping (clipped sequences present in SEQ)
        #H   hard clipping (clipped sequences NOT present in SEQ)
        #P   padding (silent deletion from padded reference)
        #=   sequence match
        #X   sequence mismatch
        """
    pairs = []
    last_position = 0
    current_position = 0
    for c in cigar:
        if not c.isdigit():
            pairs.append((c,int(cigar[last_position:current_position])))
            last_position = current_position+1
        current_position += 1

    for i in pairs:
        print(i)
    # return seq

def align_read(query, reference, matrix, gap_open=20, gap_extension=1):

    result = parasail.sw_trace_striped_8(query, reference, gap_open, gap_extension, matrix)
    traceback = result.get_traceback('|', '.', ' ')
    query_start = result.cigar.beg_query
    reference_start = result.cigar.beg_ref
    cigar = result.cigar.decode.decode("UTF-8")
    print("cigar", cigar)

    result = parasail.sw_stats_striped_8(query, reference, gap_open, gap_extension, matrix)
    
    return {
        "query_start": query_start,
        "reference_start": reference_start,
        "identity": result.matches / result.len_ref,
        "cigar": cigar,
        "ref": traceback.ref,
        "comp": traceback.comp,
        "query": traceback.query
    }

def process_file(reads,reference):
    nuc_matrix = parasail.matrix_create("ACGT", 2, -1)
    for record in SeqIO.parse(reads, "fastq"):
        print("the record id is",record.id)
        print("the seq is",record.seq)
        stats = align_read(str(record.seq), reference, nuc_matrix, gap_open=10, gap_extension=1)
        print('\n',stats["query"],'\n', stats["comp"],'\n', stats["ref"], stats["query_start"], stats["reference_start"])
        print(stats["ref"][61])
        get_site(64, stats)

        # pairs = parse_cigar(stats["cigar"], reference)

def get_site(cpg_index, stats):
    adjusted_index = cpg_index - 1 - stats["reference_start"]
    current_index = 0

    for i in range(len(stats["ref"])):
        if adjusted_index == current_index:
            print("CPG site:")
            print(stats["ref"][i],stats["query"][i])

        if stats["ref"][i] != '-':
            current_index +=1




if __name__ == '__main__':

    args = parse_args()
    ref_seq = ""
    for record in SeqIO.parse(args.reference, "fasta"):
        ref_seq = str(record.seq)
    process_file(str(args.reads), ref_seq)