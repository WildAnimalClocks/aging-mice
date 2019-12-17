from Bio import SeqIO
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Gathering and filtering files by read length.')

    parser.add_argument("--min_length", action="store", type=int, dest="min_length")
    parser.add_argument("--max_length", action="store", type=int, dest="max_length")

    parser.add_argument("--path_to_fastq", action="store", type=str, dest="path_to_fastq")
    parser.add_argument("--output_file", action="store", type=str, dest="output_file")

    return parser.parse_args()

def count_files(directory):
    file_count = 0

    for r,d,files in os.walk(directory):
        
        for f in files:
            if f.endswith(".fastq") or f.endswith(".fq"):
                file_count +=1
    return file_count

if __name__ == '__main__':

    args = parse_args()

    with open(str(args.output_file), "w") as fw:

        file_count = count_files(str(args.path_to_fastq))

        if file_count >=1 :
            print(f"Found {file_count} fastq files.\nGathering now.")

            duplicates = set()
            unique_reads = 0
            total_reads = 0 
            
            records = []

            for r,d,files in os.walk(str(args.path_to_fastq)):

                current_dir = r.split('/')[-1]
                fastq_files = [i for i in files if i.endswith(".fastq") or i.endswith("fq")]

                print(f"\tProcessing {len(fastq_files)} files in {current_dir}")

                for f in fastq_files:
                    for record in SeqIO.parse(r + '/' + f, "fastq"):
                        if len(record) > int(args.max_length):
                            continue
                        if len(record) < int(args.min_length):
                            continue

                        total_reads += 1
                        if record.id not in duplicates:
                            records.append(record)

                            duplicates.add(record.id)

                            unique_reads += 1

            SeqIO.write(records, fw, "fastq")
            print(f"Write {unique_reads} reads to {args.output_file}.")
            dups = total_reads - unique_reads
            print(f"Found {dups} duplicate reads.")

        else:
            print(f"Found {file_count} fastq files in this directory.\nExiting.")


