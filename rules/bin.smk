# rule map_to_assign_gene:
#     input:
#         reads = config["output_path"]+ "/demultiplexed_reads/{barcode}.fastq",
#         ref = config["references_file"]
#     output:
#         config["output_path"]+ "/mapping_information/{barcode}.paf"
#     shell:
#         "minimap2 -x map-ont --secondary=no --paf-no-hit {input.ref} {input.reads} > {output}"

# rule parse_paf_files:
#     input:
#         ref = config["references_file"],
#         reads = config["output_path"]+ "/demultiplexed_reads/{barcode}.fastq",
#         paf = rules.map_to_assign_gene.output
#     params:
#         path_to_script = workflow.current_basedir,
#         output_dir = config["output_path"] + "/binned/{barcode}",
#         barcode = "{barcode}"
#     output:
#         expand(config["output_path"] + "/binned/{{barcode}}/{gene}.fastq", gene=config["genes"]),
#         temp_report = temp(config["output_path"] + "/binned/{barcode}/report.csv")
#     shell:
#         """
#         python {params.path_to_script}/../scripts/parse_paf.py \
#         --paf_file {input.paf} \
#         --reads {input.reads} \
#         --output_dir {params.output_dir} \
#         --report {output.temp_report} \
#         --barcode {params.barcode} \
#         --references {input.ref}
#         """

rule collate_reports:
    input:
        expand(config["output_path"] + "/binned/{barcode}/report.csv", barcode=config["barcodes"])
    output:
        config["output_path"] + "/binned/report.csv"
    run:
        with open(str(output),"w") as fw:
            fw.write("barcode,reference,mapping_count\n")
            for input_file in input:
                with open(input_file, "r") as f:
                    for l in f:
                        l = l.rstrip()
                        fw.write(l + '\n')

rule makeblastdb:
    input:
        config["references_file"]
    output:
        config["references_file"] + ".nhr"
    shell:
        "makeblastdb -in {input} -dbtype nucl"

rule fastq_to_fasta:
    input:
        config["output_path"]+ "/demultiplexed_reads/{barcode}.fastq"
    output:
        config["output_path"]+ "/blast/{barcode}.fasta"
    shell:
        "seqtk seq -A {input} > {output}"

rule blastn:
    input:
        db=config["references_file"],
        db_hidden= rules.makeblastdb.output,
        reads=rules.fastq_to_fasta.output
    output:
        config["output_path"] + "/blast/{barcode}_blast.csv"
    shell:
        "blastn -task blastn -db {input.db} "
        "-query {input.reads} -out {output} "
        "-num_threads 1 -outfmt 10"
            
rule assign_genes:
    input:
        blast=rules.blastn.output,
        reads=config["output_path"]+ "/demultiplexed_reads/{barcode}.fastq",
        references=config["references_file"]
    params:
        outdir=config["output_path"]+ "/binned/{barcode}",
        sample="{barcode}"
    output:
        summary=temp(config["output_path"] + "/binned/{barcode}/report.csv"),
        reads=expand(config["output_path"] + "/binned/{{barcode}}/{gene}.fastq", gene=config["genes"])
    shell:
        """
        python scripts/bin.py \
        --blast_file {input.blast} \
        --references {input.references} \
        --reads {input.reads} \
        --output_dir {params.outdir} \
        --report {output.summary} \
        --barcode {params.sample} 
        """
