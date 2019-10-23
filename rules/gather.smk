rule gather:
    input:
    params:
        run_name= config["run_name"],
        path_to_fastq= config["path_to_fastq"],
        min_length= config["min_length"],
        max_length= config["max_length"],
        outdir = config["output_directory"]
    output:
        reads="{outdir}/{run_name}.fastq"
    shell:
        """
        python scripts/gather_filter.py 
        --min_length {params.min_length} 
        --max_length {params.max_length} 
        --run_name {params.run_name} 
        --path_to_fastq {params.path_to_fastq} 
        --output_directory {params.outdir}
        """
