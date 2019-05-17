rule nanopolish_index:
    input:
        summary="pipeline_output/"+ run_name+"_sequencing_summary.txt",
        reads="pipeline_output/" + run_name+"_all.fastq"
    params:
        path_to_fast5= lambda wildcards : config["path_to_fast5"]
    output:
        fai= "pipeline_output/"+run_name+"_all.fastq.index.fai",
        gzi= "pipeline_output/"+run_name+"_all.fastq.index.gzi",
        readdb="pipeline_output/"+ run_name+"_all.fastq.index.readdb",
        index="pipeline_output/"+run_name+"_all.fastq.index"
    shell:
        "nanopolish index -d {params.path_to_fast5} {input.reads}"




