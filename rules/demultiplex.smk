rule demultiplex_porechop:
    input:
    params:
        gathered_file = config["output_path"] + "/" + config["run_name"]+ ".fastq",
        outdir = config["output_path"]+"/demultiplexed_reads/",

        discard_middle=discard_middle,
        split_reads=split_reads,
        discard_unassigned=discard_unassigned,

        require_two_barcodes=require_two_barcodes,
        barcode_option = barcode_set,
        limit_barcodes_to = limit_barcodes_to,

        threshold = "--barcode_threshold " + str(config["barcode_threshold"]),
        diff = "--barcode_diff " + str(config["barcode_diff"]),

    threads:
        2
    output:
        expand(config["output_path"]+ "/demultiplexed_reads/{barcode}.fastq", barcode = config["barcodes"])
    shell:
        """
        porechop \
        --verbosity 2 \
        -i {params.gathered_file} \
        -b {params.outdir} \
        --threads 2 \
        --barcode_labels \
        --untrimmed \
        {params.threshold} \
        {params.diff}\
        {params.limit_barcodes_to}\
        {params.require_two_barcodes}\
        {params.discard_middle}\
        {params.split_reads} \
        {params.discard_unassigned}\
        {params.barcode_option}
        """

# To do:
# add in check for when no reads are present for a particular barcode. 
# It currently would break the snakemake and you'd have to edit the config and remove that barcode.
# Fix may be to run the demultiplexing, touch the output if not there, 
# check if empty file, remove file, edit config_dict with persistent dict. 
