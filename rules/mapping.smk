# will have a sorting, qc and variant calling setup here
# nanopolish? 
# rule minimap2:
#     input:
#         fastq="pipeline_output/binned/{barcode}/reads/{gene}.fastq",
#         ref="pipeline_output/{barcode}/nanopolish_ref/{gene}.fasta"
#     output:
#         "pipeline_output/binned/{barcode}/mapped/{gene}.sam"
#     threads: 4
#     shell:
#         "minimap2 -ax map-ont --sam-hit-only --secondary=no {input.ref} {input.fastq} > {output}"

# rule samtools_sort:
#     input:
#         "pipeline_output/binned/{barcode}/mapped/{gene}.sam"
#     output:
#         "pipeline_output/binned/{barcode}/sorted_reads/{gene}.bam"
#     shell:
#         "samtools sort -T pipeline_output/binned/{wildcards.barcode}/sorted_reads/{wildcards.gene} "
#         "-O bam > {output}"

# rule samtools_index:
#     input:
#         "pipeline_output/binned/{barcode}/sorted_reads/{gene}.bam"
#     output:
#         "pipeline_output/binned/{barcode}/sorted_reads/{gene}.bam.bai"
#     shell:
#         "samtools index {input}"

# rule filter_normalise: to do- external script
#     input:
#         "pipeline_output/binned/{barcode}/sorted_reads/{gene}.bam"
#     output:
#         "pipeline_output/binned/{barcode}/filtered_reads/{gene}.bam"
#     run:
#         aln_file = str(input)


