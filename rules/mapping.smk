# will have a sorting, qc and variant calling setup here
# nanopolish? 
rule minimap2:
    input:
        fastq="pipeline_output/binned/{barcode}/reads/{gene}.fastq",
        ref="pipeline_output/{barcode}/nanopolish_ref/{gene}.fasta"
    output:
        "pipeline_output/binned/{barcode}/mapped/{gene}.sam"
    threads: 4
    shell:
        "minimap2 -ax map-ont --secondary=no {input.ref} {input.fastq} > {output}"
