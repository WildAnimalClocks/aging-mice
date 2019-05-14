rule makeblastdb:
    input:
        "references/genes.fasta"
    output:
        "references/genes.fasta.nhr"
    shell:
        "makeblastdb -in {input} -dbtype nucl"

rule fastq_to_fasta:
    input:
        "demultiplexed/{barcode}.fastq"
    output:
        "demultiplexed/{barcode}.fasta"
    shell:
        "seqtk seq -A {input} > {output}"

rule blastn:
    input:
        db_hidden="references/genes.fasta.nhr",
        db="references/genes.fasta",
        reads="demultiplexed/{barcode}.fasta" 
    output:
        "blast_results/{barcode}.blast.csv"
    shell:
        "blastn -task blastn -db {input.db} "
        "-query {input.reads} -out {output} "
        "-num_threads 16 -outfmt 10"
rule bin:
    input:
        blast="blast_results/{barcode}.blast.csv",
        reads="demultiplexed/{barcode}.fastq"
    params:
        outdir="binned/{barcode}"
    output:
        summary="binned/{barcode}/binning_report.txt"
        #fastq=expand("binned/{{barcode}}/{gene}.fastq", gene=config["genes"])
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
        for record in SeqIO.parse(reads,"fastq"):
            try:
                gene= top_gene[record.id]
            except:
                gene="none"
            records[gene].append(record)

        with open(summary,"w") as f:
            for gene in records:
                outfile = outdir + '/' + gene + '.fastq'
                f.write("{}: {} records written to {}. \n".format(gene,len(records[gene]),outfile))
                with open(outfile, "w") as fw:
                    SeqIO.write(records[gene], fw, "fastq")


                


