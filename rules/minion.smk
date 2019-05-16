rule artic_minion:
    input:
        read_file = "pipeline_output/binned/{barcode}/reads/{gene}.fastq",
        nano_read_file = "pipeline_output/"+run_name + "_all.fastq",
        index="pipeline_output/"+run_name+"_all.fastq.index"
    params:
        primer_scheme = lambda wildcards : config["primer_scheme_dir"],
        primer_version = "V_{barcode}_{gene}"
    threads:
        8
    output:
        "{barcode}_{gene}.primertrimmed.sorted.bam",
        "{barcode}_{gene}.primertrimmed.vcf",
        "{barcode}_{gene}.alignreport.txt",
        "{barcode}_{gene}.consensus.fasta"
    shell:
        "artic minion --normalise 200 --threads 16 "
        "--scheme-directory {params.primer_scheme} "
        "--read-file {input.read_file} "
        "--nanopolish-read-file {input.nano_read_file} "
        "{params.primer_version} {wildcards.barcode}"

