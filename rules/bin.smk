rule makeblastdb:
    input:
        config["references_file"]
    output:
        config["references_file"] + ".nhr"
    shell:
        "makeblastdb -in {input} -dbtype nucl"

rule fastq_to_fasta:
    input:
        config["output_path"]+ "/{barcode}_bin/{barcode}.fastq"
    output:
        config["output_path"]+ "/{barcode}_bin/{barcode}.fasta"
    shell:
        "seqtk seq -A {input} > {output}"

rule blastn:
    input:
        db_hidden=config["references_file"],
        db= rules.makeblastdb.output,
        reads=config["output_path"]+ "/{barcode}_bin/{barcode}.fasta"
    output:
        config["output_path"]+ "/{barcode}_bin/blast_results/{barcode}.blast.csv"
    shell:
        "blastn -task blastn -db {input.db} "
        "-query {input.reads} -out {output} "
        "-num_threads 1 -outfmt 10"
            
rule bin:
    input:
        blast=config["output_path"]+ "/{barcode}_bin/blast_results/{barcode}.blast.csv",
        reads=config["output_path"]+ "/{barcode}_bin/{barcode}.fastq",
        references=config["references_file"],
        primers=config["primers_file"]
    params:
        outdir=config["output_path"]+ "/{barcode}_bin/reads/",
        sample="{barcode}"
    output:
        summary=config["output_path"]+ "/{barcode}_bin/binning_report.txt",
        ref=expand(config["output_path"] + "/{{barcode}}_bin/{gene}.reference.fasta", gene=config["genes"]),
        bed=expand(config["output_path"] + "/{{barcode}}_bin/{gene}.primers.bed", gene=config["genes"]),
        reads=expand(config["output_path"]+ "/{{barcode}}_bin/reads/{gene}.fastq", gene=config["genes"])
    shell:
        """
        python scripts/bin.py 
        --blast_file {input.blast} 
        --reference_file {input.references} 
        --reads {input.reads} 
        --output_dir {params.outdir} 
        --summary {output.summary} 
        --sample {params.sample} 
        --primers {input.primers}
        """
