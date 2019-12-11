rule chorepop_align:
    input:
        amp_seqs = "references/genes.modified.converted.fasta",
        reads = config["output_path"] + "/binned/{{barcode}}/{gene}.fastq"
    output:
        

# rule align:
#     input:
#         cns="pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.consensus.fasta",
#         refs="references/mod_and_unmod/{gene}.fasta"
#     output:
#         "pipeline_output/mod_counting/{barcode}_bin/{barcode}_{gene}.aln.fasta"
#     shell:
#         "cat {input.cns} {input.refs} > pipeline_output/mod_counting/{barcode}_bin/cns_and_ref_{wildcards.gene}.fasta && "
#         "mafft pipeline_output/mod_counting/{barcode}_bin/cns_and_ref_{wildcards.gene}.fasta > {output}"

rule count:
    input:
        vcf=expand("pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.vcf",barcode=config["barcodes"], gene=config["genes"]),
        cpg="references/cpg_sites.csv",
        age="references/ages.csv"
    output:
        cpg="pipeline_output/cpg_report.csv",
        cpg_wide="pipeline_output/cpg_per_position.csv"
    shell:
        "python scripts/variant_counts.py --cpg_info {input.cpg} --out_file {output.cpg} --ages {input.age} --cpg_wide {output.cpg_wide}"
