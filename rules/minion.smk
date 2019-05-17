rule artic_minion:
    input:
        read_file = "pipeline_output/binned/{barcode}_bin/reads/{gene}.fastq",
        nano_read_file = "pipeline_output/"+run_name + "_all.fastq",
        index="pipeline_output/"+run_name+"_all.fastq.index",
        fasta="pipeline_output/binned/{barcode}_bin/primer-schemes/minion/V_{gene}/minion.reference.fasta",
        bed="pipeline_output/binned/{barcode}_bin/primer-schemes/minion/V_{gene}/minion.scheme.bed"
    params:
        primer_scheme = "pipeline_output/binned/{barcode}_bin/primer-schemes",
        primer_version = "minion/V_{gene}",
        sample = "{barcode}_{gene}"
    threads:
        8
    output:
        "{barcode}_{gene}.alignreport.er",
        "{barcode}_{gene}.alignreport.txt",
        "{barcode}_{gene}.consensus.fasta",
        "{barcode}_{gene}.minion.log.txt",
        "{barcode}_{gene}.primertrimmed.sorted.bam",
        "{barcode}_{gene}.primertrimmed.sorted.bam.bai",
        "{barcode}_{gene}.primertrimmed.vcf",
        "{barcode}_{gene}.sorted.bam",
        "{barcode}_{gene}.sorted.bam.bai",
        "{barcode}_{gene}.trimmed.sorted.bam",
        "{barcode}_{gene}.trimmed.sorted.bam.bai",
        "{barcode}_{gene}.variants.tab",
        "{barcode}_{gene}.vcf"
    shell:
        "artic minion --normalise 200 --threads 1 "
        "--scheme-directory {params.primer_scheme} "
        "--read-file {input.read_file} "
        "--nanopolish-read-file {input.nano_read_file} "
        "{params.primer_version} {params.sample}"
        
rule organise_minion_output:
    input:
        "{barcode}_{gene}.alignreport.er",
        "{barcode}_{gene}.alignreport.txt",
        "{barcode}_{gene}.consensus.fasta",
        "{barcode}_{gene}.minion.log.txt",
        "{barcode}_{gene}.primertrimmed.sorted.bam",
        "{barcode}_{gene}.primertrimmed.sorted.bam.bai",
        "{barcode}_{gene}.primertrimmed.vcf",
        "{barcode}_{gene}.sorted.bam",
        "{barcode}_{gene}.sorted.bam.bai",
        "{barcode}_{gene}.trimmed.sorted.bam",
        "{barcode}_{gene}.trimmed.sorted.bam.bai",
        "{barcode}_{gene}.variants.tab",
        "{barcode}_{gene}.vcf"
    params:
        output_dir="pipeline_output/minion_output/{barcode}"
    output:
        "pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.alignreport.er",
        "pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.alignreport.txt",
        "pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.consensus.fasta",
        "pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.minion.log.txt",
        "pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.primertrimmed.sorted.bam",
        "pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.primertrimmed.sorted.bam.bai",
        "pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.primertrimmed.vcf",
        "pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.sorted.bam",
        "pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.sorted.bam.bai",
        "pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.trimmed.sorted.bam",
        "pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.trimmed.sorted.bam.bai",
        "pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.variants.tab",
        "pipeline_output/minion_output/{barcode}_bin/{barcode}_{gene}.vcf"
    shell:
        "mv {input} {params.output_dir}"
