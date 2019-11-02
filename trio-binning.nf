#!/usr/bin/env nextflow

// pachon is maternal; surface is paternal

params.maternal_short_reads = 'reads/short_reads/AGTTCGTC-*.fastq.gz'
params.paternal_short_reads = 'reads/short_reads/GTCATCGA-*.fastq.gz'
params.input_reads = 'reads/long_reads/*.subreads.bam'
params.genomeSize = '1.6g'

maternal_short_reads = Channel.fromPath(params.maternal_short_reads)
paternal_short_reads = Channel.fromPath(params.paternal_short_reads)
input_reads = Channel
    .fromPath(params.input_reads)
    .map { f -> tuple(f.baseName, f) }

process find_unique_kmers {
    cpus 16
    publishDir 'unique_kmers'

    input:
    file maternal from maternal_short_reads.collect()
    file paternal from paternal_short_reads.collect()

    output:
    file 'hapA_only_kmers.txt' into maternal_kmers
    file 'hapB_only_kmers.txt' into paternal_kmers

    """
    find-unique-kmers -k 20 -p 16 ${maternal.join(',')} ${paternal.join(',')}
    """
}

process classify_long_reads {
    cpus 2
    clusterOptions '--qos long'
    time '7d'
    publishDir 'classified_reads'

    input:
    file 'maternal_kmers.txt' from maternal_kmers
    file 'paternal_kmers.txt' from paternal_kmers
    set id, file(reads_to_classify) from input_reads

    output:
    file "maternal.${id}.fq.gz" into maternal_reads
    file "paternal.${id}.fq.gz" into paternal_reads
    file "unclassified.${id}.fq.gz" into unclassified_reads

    """
    classify_by_kmers -c -i $reads_to_classify \
        -a maternal_kmers.txt -b paternal_kmers.txt \
        -A maternal.${id} -B paternal.${id} -U unclassified.${id}
    """
}

maternal_reads.map { tuple('maternal', it) }.set { maternal_tuples }
paternal_reads.map { tuple('paternal', it) }.set { paternal_tuples }

maternal_tuples.concat(paternal_tuples).groupTuple().set { both_haplotypes }

process assemble {
    clusterOptions '--qos long'
    time '7d'
    memory '12g'
    publishDir 'contigs'

    input:
    set haplotype, file("${haplotype}.*.fq.gz") from both_haplotypes

    output:
    file "${haplotype}/${haplotype}.contigs.fasta" into contigs

    """
    /storage/htc/warrenlab/esrbhb/software/canu-1.8/Linux-amd64/bin/canu \
        -p ${haplotype} -d ${haplotype} \
        genomeSize=${params.genomeSize} \
        gridOptions="-p Lewis,BioCompute -t 2-00:00:00 --account warrenlab" \
        -pacbio-raw ${haplotype}.*.fq.gz
    """
}
