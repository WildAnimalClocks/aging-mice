# will have a sorting, qc and variant calling setup here
# nanopolish? 

rule minimap2:
    input:
        fastq="binned/{barcode}/{gene}.fastq",
        ref="references/genes/{gene}.fasta"
    output:
        "binned/{barcode}/mapped/{gene}.sam"
    threads: 4
    shell:
        "minimap2 -ax map-ont --secondary=no {input.ref} {input.fastq} > {output}"  