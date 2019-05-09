rule demultiplex_qcat:
    input:
        reads=expand(run_name + "_all.fastq")
    output:
        fastq=expand("demultiplexed/{barcode}.fastq",barcode=config["barcodes"]),
        report="data/demultiplex_report.txt"
    threads: 15
    shell:
        "qcat -f {input.reads} -b data/demultiplexed -t 16 --min-score 80 > {output.report}"

# rule demultiplex_porechop:
#     input:
#         reads=expand(run_name + "_all.fastq")
#     output:
#         fastq=expand("demultiplexed/{barcode}.fastq",barcode=config["barcodes"]),
#         report="data/demultiplex_report.txt"
#     threads: 15
#     shell:
#         "porechop -i {input.reads} --verbosity 2 --untrimmed --discard_middle "
#         "--native_barcodes --barcode_threshold 80 "
#         "--threads 16 --check_reads 10000 --barcode_diff 5 "
#         "-b data/demultiplexed > {output.report}"
