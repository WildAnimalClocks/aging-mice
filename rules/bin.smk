rule makeblastdb:
    input:
        "references/genes.all.fasta"
    output:
        "references/genes.all.fasta.nhr"
    shell:
        "makeblastdb -in {input} -dbtype nucl"

rule fastq_to_fasta:
    input:
        "pipeline_output/demultiplexed/{barcode}.fastq"
    output:
        "pipeline_output/demultiplexed/{barcode}.fasta"
    shell:
        "seqtk seq -A {input} > {output}"

rule blastn:
    input:
        db_hidden="references/genes.all.fasta.nhr",
        db="references/genes.all.fasta",
        reads="pipeline_output/demultiplexed/{barcode}.fasta" 
    output:
        "pipeline_output/blast_results/{barcode}.blast.csv"
    shell:
        "blastn -task blastn -db {input.db} "
        "-query {input.reads} -out {output} "
        "-num_threads 1 -outfmt 10"
            
rule bin:
    input:
        blast="pipeline_output/blast_results/{barcode}.blast.csv",
        reads="pipeline_output/demultiplexed/{barcode}.fastq",
        references="references/genes.all.fasta"
    params:
        outdir="pipeline_output/binned/{barcode}_bin/reads/",
        sample="{barcode}"
    output:
        summary="pipeline_output/binned/{barcode}_bin/binning_report.txt",
        ref=expand("pipeline_output/binned/{{barcode}}_bin/primer-schemes/minion/V_{gene}/minion.reference.fasta", gene=config["genes"]),
        bed=expand("pipeline_output/binned/{{barcode}}_bin/primer-schemes/minion/V_{gene}/minion.scheme.bed", gene=config["genes"]),
        reads=expand("pipeline_output/binned/{{barcode}}_bin/reads/{gene}.fastq", gene=config["genes"])
    shell:
        "python scripts/bin.py --blast_file {input.blast} "
        "--reference_file {input.references} --reads {input.reads} "
        "--output_dir {params.outdir} --summary {output.summary} "
        "--sample {params.sample}"
