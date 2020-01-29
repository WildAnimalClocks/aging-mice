import argparse
from Bio import SeqIO
import csv
from collections import defaultdict
from collections import Counter

def parse_args():
    parser = argparse.ArgumentParser(description='Parse barcode info from reads and bin.')

    parser.add_argument("--output_prefix", action="store", type=str, dest="output_prefix")
    parser.add_argument("--annotated_reads", action="store", type=str, dest="annotated_reads")
    parser.add_argument("--barcode_csv", action="store", type=str, dest="barcode_csv")

    parser.add_argument("--barcode_threshold", default=80, action="store", type=float, dest="barcode_threshold")

    return parser.parse_args()

def parse_reference_file(references):
    #returns a dict of dicts containing reference header information
    #key is seq id i.e. first field of the header string
    ref_info = defaultdict(dict)
    for record in SeqIO.parse(references,"fasta"):
        tokens = record.description.split(' ')
        for i in tokens:
            try:
                info = i.split('=')
                ref_info[record.id][info[0]]=info[1]
                # ref_info['*'][info[0]]='NA'
            except:
                pass
    return ref_info

def parse_read_header(header):
    #returns a dict of {key:value} pairs containing all 
    #" key=value" strings present on the read header

    tokens= header.split(' ')
    header_info = {}
    for i in tokens:
        try:
            info = i.split('=')
            header_info[info[0]]=info[1]
        except:
            pass
    return header_info

def get_header_dict(reads):
    #This function parses the fastq file and returns a dictionary
    #with read name as the key and barcode information as the value
    # i.e. barcode_dict[read_name]=barcode

    header_dict = {}
    for record in SeqIO.parse(str(reads),"fastq"):
        header = parse_read_header(str(record.description))
        try:
            start_name = header["start_barcode"]
            start_identity = header["start_id"]
            end_name = header["end_barcode"]
            end_identity = header["end_id"]
            start_time = header["start_time"]
        except:
            start_name = header["start_barcode"]
            start_identity = header["start_id"]
            end_name = header["end_barcode"]
            end_identity = header["end_id"]
            start_time = header["start_time"]
            
            start_time = header["start_time"]

        header_dict[record.id]=(barcode, start_time)
        
    return header_dict

def check_identity_threshold(mapping, min_identity):
    
    if float(min_identity)<1:
        min_id = float(min_identity)
    else:
        min_id = float(min_identity)/ 100
        
    if mapping["identity"] >= min_id:
        return True
    else:
        return False

def parse_line(line, header_dict):
    values = {}
    tokens = line.rstrip('\n').split('\t')
    values["read_name"], values["read_len"] = tokens[:2]
    if values["read_name"] in header_dict:
        values["barcode"], values["start_time"] = header_dict[values["read_name"]] #if porechop didn't discard the read
    else:
        values["barcode"], values["start_time"] = "none", "?" #don't have info on time or barcode
    values["query_start"] = tokens[2]
    values["query_end"] = tokens[3]
    values["ref_hit"], values["ref_len"], values["coord_start"], values["coord_end"], values["matches"], values["aln_block_len"] = tokens[5:11]
    if values["ref_hit"] != "*":
        mismatches, identity = calculate_genetic_identity(tokens[-1])
        values["mismatches"] = mismatches
        values["identity"] = identity
    else:
        values["mismatches"] = 0
        values["identity"]= 0

    return values


def write_mapping(report, mapping, reference_options, reference_info, counts, min_identity):
    if mapping["ref_hit"] == '*' or mapping["ref_hit"] == '?':
        # '*' means no mapping, '?' ambiguous mapping (i.e., multiple primary mappings)
        mapping['coord_start'], mapping['coord_end'] = 0, 0
        if (mapping["ref_hit"] == '*'):
            counts["unmapped"] += 1
        else:
            counts["ambiguous"] += 1

        if reference_options != None:
                mapping["ref_opts"] = []
                for k in reference_options:
                    mapping["ref_opts"].append(mapping["ref_hit"])
    else:
        if reference_options != None:
            mapping["ref_opts"] = []
            for k in reference_options:
                if len(reference_options[k]) == 1:
                    mapping["ref_opts"].append(reference_info[mapping["ref_hit"]][k])
                else:
                    overlap_list = []
                    for i in reference_options[k]:
                        if len(i) == 3:
                            sub_k, opt_start, opt_end = i
                            overlap, length = check_overlap((opt_start, opt_end),(int(mapping["coord_start"]), int(mapping["coord_end"])))
                            if overlap:
                                overlap_list.append((reference_info[mapping["ref_hit"]][sub_k], length))
                    if overlap_list:
                        best = sorted(overlap_list, key = lambda x : x[1], reverse=True)[0]
                        mapping["ref_opts"].append(best[0])
                    else:
                        mapping["ref_opts"].append("NA")

    

    if check_identity_threshold(mapping, min_identity):

        counts["total"] += 1

        mapping_length = int(mapping['matches']) + int(mapping['mismatches'])
        report.write(f"{mapping['read_name']},{mapping['read_len']},{mapping['start_time']},"
                    f"{mapping['barcode']},{mapping['ref_hit']},{mapping['ref_len']},"
                    f"{mapping['coord_start']},{mapping['coord_end']},{mapping['matches']},{mapping_length}")
        if 'ref_opts' in mapping:
            report.write(f",{','.join(mapping['ref_opts'])}\n")
        else:
            report.write("\n")
    else:
        counts["unmapped"] +=1
        report.write(f"{mapping['read_name']},{mapping['read_len']},{mapping['start_time']},"
                    f"{mapping['barcode']},*,0,0,0,0,0")
        if 'ref_opts' in mapping:
            ref_opt_list = ['*' for i in mapping['ref_opts']]
            report.write(f",{','.join(ref_opt_list)}\n")
        else:
            report.write("\n")

def parse_paf(paf, report, header_dict, reference_options, reference_info,min_identity):
    #This function parses the input paf file 
    #and outputs a csv report containing information relevant for RAMPART and barcode information
    # read_name,read_len,start_time,barcode,best_reference,start_coords,end_coords,ref_len,matches,aln_block_len,ref_option1,ref_option2
    counts = {
        "unmapped": 0,
        "ambiguous": 0,
        "total": 0
    }

    with open(str(paf),"r") as f:
        last_mapping = None
        for line in f:

            mapping = parse_line(line, header_dict)

            if last_mapping:
                if mapping["read_name"] == last_mapping["read_name"]:
                    # this is another mapping for the same read so set the original one to ambiguous. Don't
                    # set last_mapping in case there is another mapping with the same read name.
                    last_mapping['ref_hit'] = '?'
                else:
                    write_mapping(report, last_mapping, reference_options, reference_info, counts,min_identity)
                    last_mapping = mapping
            else:
                last_mapping = mapping

        # write the last last_mapping
        write_mapping(report, last_mapping, reference_options, reference_info, counts,min_identity)

    try:
        prop_unmapped = counts["unmapped"] / counts["total"]
        print("Proportion unmapped is {}".format(prop_unmapped))
        if prop_unmapped >0.95:
            print("\nWarning: Very few reads have mapped (less than 5%).\n")
    except:
        print("Probably can't find the records.") #division of zero the error

if __name__ == '__main__':

    args = parse_args()

    with open(str(args.report), "w") as csv_report:
        
        if args.reference_options:
            reference_options,ref_option_header = parse_reference_options(args.reference_options)
            reference_info = parse_reference_file(args.references)
        else:
            reference_options,ref_option_header = None,''
            reference_info = None

        header_dict = get_header_dict(args.reads)

        csv_report.write(f"read_name,read_len,start_time,barcode,best_reference,ref_len,start_coords,end_coords,num_matches,mapping_len{ref_option_header}\n")
        parse_paf(args.paf_file, csv_report, header_dict, reference_options, reference_info,args.min_identity)
