rule gather:
    input:
    params:
        run_name= config["run_name"],
        path_to_fastq= config["input_path"],
        min_length= config["min_length"],
        max_length= config["max_length"],
        outdir = config["output_path"],
        path_to_script = workflow.current_basedir
    output:
        reads="{outdir}/{run_name}.fastq"
    shell:
        """
        python {params.path_to_script}/../scripts/gather_filter.py \
        --min_length {params.min_length} \
        --max_length {params.max_length} \
        --path_to_fastq {params.path_to_fastq} \
        --output_file {output.reads} 
        """
