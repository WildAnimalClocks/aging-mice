rule makeblastdb:
    input:
        "references/genes.fasta"
    output:
        "references/genes.fasta"
    shell:
        "makeblastdb -in {input} -dbtype nucl"

rule fastq_to_fasta:
    input:
        "demultiplexed/{barcode}.fastq"
    output:
        "demultiplexed/{barcode}.fasta"
    shell:
        "seqtk seq -A {input} > {output}"

rule blastn:
    input:
        db_hidden="references/genes.fasta.nhr",
        db="references/genes.fasta",
        reads="demultiplexed/{barcode}.fasta" 
    output:
        "blast_results/{barcode}.blast.txt"
    shell:
        "blastn -task blastn -db {input.db} "
        "-query {input.reads} -out {output} "
        "-num_threads 4 -outfmt 6"
rule bin:
    input:
        blast="blast_results/{barcode}.blast.txt",
        reads="demultiplexed/{barcode}.fastq"
    output:
        summary="blast_results/parsed_blast_info.csv",
        fastq=expand("binned/{barcode}/{gene}.fastq",barcode=config["barcodes"], gene=config["genes"])
    script:


