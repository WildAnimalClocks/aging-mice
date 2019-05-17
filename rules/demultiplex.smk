rule demultiplex_qcat:
    input:
        reads="pipeline_output/"+run_name + "_all.fastq"
    output:
        fastq=expand("pipeline_output/demultiplexed/{barcode}.fastq",barcode=config["barcodes"]),
        report="pipeline_output/demultiplexed/demultiplex_report.txt"
    threads: 8
    shell:
        "qcat -f {input.reads} -b pipeline_output/demultiplexed -t 8 -q 80 > {output.report}"

