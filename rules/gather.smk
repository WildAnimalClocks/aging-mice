rule gather:
    input:
    params:
        run_name=lambda wildcards : {config["run_name"]},
        path_to_fastq= lambda wildcards : {config["path_to_fastq"]},
        min_length= lambda wildcards : {config["min_length"]},
        max_length= lambda wildcards : {config["max_length"]}
    output:
        summary="pipeline_output/{run_name}_sequencing_summary.txt",
        reads="pipeline_output/{run_name}_all.fastq"
    shell:
        "python scripts/gather_filter.py --min-length {params.min_length} "
        "--max-length {params.max_length} --run-name {params.run_name} "
        "--path-to-fastq {params.path_to_fastq}"
