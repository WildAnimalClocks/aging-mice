rule demultiplex_qcat:
    input:
        reads=run_name + "_all.fastq"
    output:
        fastq=expand("demultiplexed/{barcode}.fastq",barcode=config["barcodes"]),
        report="demultiplexed/demultiplex_report.txt"
    threads: 16
    shell:
        "qcat -f {input.reads} -b demultiplexed -t 16 -q 80 > {output.report}"

# rule demultiplex_porechop:
#     input:
#         reads=expand(run_name + "_all.fastq")
#     output:
#         fastq=expand("demultiplexed/{barcode}.fastq",barcode=config["barcodes"]),
#         report="demultiplexed/demultiplex_report.txt"
#     threads: 15
#     shell:
#         "porechop -i {input.reads} --verbosity 2 --untrimmed --discard_middle "
#         "--native_barcodes --barcode_threshold 80 "
#         "--threads 16 --check_reads 10000 --barcode_diff 5 "
#         "-b data/demultiplexed > {output.report}"
