rule makeblastdb:
    input:
        "references/genes.all.fasta"
    output:
        "references/genes.all.fasta.nhr"
    shell:
        "makeblastdb -in {input} -dbtype nucl"

rule fastq_to_fasta:
    input:
        "pipeline_output/demultiplexed/{barcode}.fastq"
    output:
        "pipeline_output/demultiplexed/{barcode}.fasta"
    shell:
        "seqtk seq -A {input} > {output}"

rule blastn:
    input:
        db_hidden="references/genes.all.fasta.nhr",
        db="references/genes.all.fasta",
        reads="pipeline_output/demultiplexed/{barcode}.fasta" 
    output:
        "pipeline_output/blast_results/{barcode}.blast.csv"
    shell:
        "blastn -task blastn -db {input.db} "
        "-query {input.reads} -out {output} "
        "-num_threads 16 -outfmt 10"
            
rule bin:
    input:
        blast="pipeline_output/blast_results/{barcode}.blast.csv",
        reads="pipeline_output/demultiplexed/{barcode}.fastq",
        references="references/genes.all.fasta"
    params:
        outdir="pipeline_output/binned/{barcode}/reads"
    output:
        summary="pipeline_output/binned/{barcode}/binning_report.txt",
        ref=expand("pipeline_output/{{barcode}}/nanopolish_ref/{gene}.fasta", gene=config["genes"]),
        fastq=expand("pipeline_output/binned/{{barcode}}/reads/{gene}.fastq", gene=config["genes"])
    run:
        blast_file = str(input.blast)
        reads = str(input.reads)
        outdir = str(params.outdir)
        summary = str(output.summary)

        hits = collections.defaultdict(list)
        with open(blast_file, "r") as f:
            for l in f:
                tokens= l.rstrip('\n').split(',')
                read = tokens[0]
                hit = tokens[1]
                score= tokens[-1] #check these
                hits[read].append((hit,score))

        top_gene = {}
        for i in hits:
            top_hit = sorted(hits[i], key = lambda x : float(x[1]), reverse = True)[0]
            top_gene[i]=top_hit[0]
        hits.clear()

        records = collections.defaultdict(list)
        mod_counter = collections.Counter()
        for record in SeqIO.parse(reads,"fastq"):
            try:
                gene,mod= top_gene[record.id].split('_')
                mod_counter[top_gene[record.id]]+=1
            except:
                gene="none"
                mod_counter["none"]+=1
            records[gene].append(record)

        seq_dict = {}
        for record in SeqIO.parse(str(input.references),"fasta"):
            seq_dict[record.id]=record.seq.upper()

        with open(summary,"w") as f:
            for gene in records:
                unmodified_count = mod_counter[gene+"_unmodified"]
                modified_count = mod_counter[gene+"_modified"]
                if unmodified_count >= modified_count:
                    #use the unmodified gene as reference for mapping
                    with open("pipeline_output/nanopolish_ref/"+gene+".fasta","w") as fnano:
                        fnano.write(">{}_unmodified\n{}\n".format(gene,seq_dict[gene+"_unmodified"]))
                else:
                    with open("pipeline_output/nanopolish_ref/"+gene+".fasta","w") as fnano:
                        fnano.write(">{}_modified\n{}\n".format(gene,seq_dict[gene+"_modified"]))

                print("In file:{}\n For gene: {}\n\tUnmodified hit count: {}\n\tModified hit count: {}\n".format(blast_file, gene, unmodified_count, modified_count))
                out_file = outdir + '/' + gene + '.fastq'
                f.write("{}: {} records written to {}. \n".format(gene,len(records[gene]),out_file))
                with open(out_file, "w") as fw:
                    SeqIO.write(records[gene], fw, "fastq")

