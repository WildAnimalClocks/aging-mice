rule ParaMethR:
    input:
        refs = "references/genes.converted.ambiguous.fasta",
        reads = config["output_path"] + "/demultiplexed_reads/{barcode}.fastq",
        cpg = "references/cpg_sites.csv",
        matrix = "references/substitution_matrix.txt",
    params:
        path_to_script = workflow.current_basedir,
        sample = "{barcode}"
    output:
        report = config["output_path"] +"/processed_reads/counts/{barcode}.cpg_counts.csv",
        counts = config["output_path"] +"/processed_reads/counts/{barcode}.gene_counts.csv"
    shell:
        """
        python {params.path_to_script}/../scripts/paramethr.py \
            --reads {input.reads} \
            --reference {input.refs} \
            --cpg_csv {input.cpg} \
            --substitution_matrix {input.matrix} \
            --sample {params.sample} \
            --report {output.report} \
            --counts {output.counts}
        """

rule gather_reports:
    input:
        expand(config["output_path"] +"/processed_reads/counts/{barcode}.gene_counts.csv", barcode = config["barcodes"])
    output:
        config["output_path"] +"/reports/cpg_counts.csv"
    shell:
        "cat {input} > {output}"

rule gather_counts:
    input:
        expand(config["output_path"] +"/processed_reads/counts/{barcode}.cpg_counts.csv", barcode = config["barcodes"])
    output:
        config["output_path"] +"/reports/gene_counts.csv"
    shell:
        "cat {input} > {output}"