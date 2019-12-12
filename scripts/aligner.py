import argparse
import os
import sys
from collections import Counter
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

def get_best_reference(query, ref_dict, matrix):
    best_reference_alignment = {
        "reference": "None",
        "identity": 0
        }
    for ref in ref_dict:
        new_reference_alignment = align_read(query, ref, ref_dict[ref], matrix)

        if best_reference_alignment["identity"] < new_reference_alignment["identity"]:
            best_reference_alignment = new_reference_alignment
        else:
            pass
    return best_reference_alignment

def align_read(query, ref_id, reference, matrix, gap_open=20, gap_extension=1):

    result = parasail.sw_trace_striped_8(query, reference, gap_open, gap_extension, matrix)
    traceback = result.get_traceback('|', '.', ' ')
    query_start = result.cigar.beg_query
    reference_start = result.cigar.beg_ref
    # cigar = result.cigar.decode.decode("UTF-8")

    result = parasail.sw_stats_striped_8(query, reference, gap_open, gap_extension, matrix)
    try:
        return {
            "reference":ref_id,
            "query_start": query_start,
            "reference_start": reference_start,
            "matches": result.matches,
            "len": result.len_ref,
            "identity": result.matches / result.len_ref,
            "ref": traceback.ref,
            "comp": traceback.comp,
            "query": traceback.query
        }
    except:
        return {
            "reference":ref_id,
            "identity":0,
            "query_start": query_start,
            "len": result.len_ref,
            "reference_start": reference_start,
            "ref": traceback.ref,
            "comp": traceback.comp,
            "query": traceback.query
        }

def process_file(reads,references):

    nuc_matrix = parasail.matrix_create("ACGT", 2, -1)
    counts = Counter()
    
    for record in SeqIO.parse(reads, "fastq"):

        # print("*****")
        # print("the record id is",record.id)
        stats = get_best_reference(str(record.seq), references, nuc_matrix)
        if stats["identity"] > 0.5:

            print('\n', stats["reference"], '\n', stats["query"],'\n', stats["comp"],'\n', stats["ref"])
            print(record.id, stats["identity"], stats["query_start"], stats["reference_start"])
            counts[stats["reference"]]+=1

        else:
            counts["no homology"]+=1
    return counts
        # print(stats["ref"][61])
        # get_site(64, stats)


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
    references = {}
    for record in SeqIO.parse(args.reference, "fasta"):
         references[record.id] = str(record.seq)

    counts = process_file(str(args.reads), references)
    for i in counts:
        print(i, '\t', counts[i])