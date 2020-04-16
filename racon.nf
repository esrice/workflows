#!/usr/bin/env nextflow

params.reference = '../contigs/jackson-trio/jackson-haplotypematernal/jackson-haplotypematernal.contigs.fasta'
params.reads = '../contigs/reads/long_reads/*.fastq.gz'
params.numChunks = 100

reference = file(params.reference)
Channel.fromPath(params.reads)
    .map { f -> tuple(f.baseName, f) }
    .set { readsToAlign }
Channel.fromPath(params.reads)
    .set { readsToCat }

process faidx {
    input:
    file 'ref.fa' from reference

    output:
    file 'ref.fa.fai' into referenceFaidx

    """
    samtools faidx ref.fa
    """
}

process align {
    cpus 16

    input:
    file 'ref.fa' from reference
    set id, file(inputReadFastq) from readsToAlign

    output:
    file "${id}.paf" into aligned

    """
    minimap2 -x ava-ont -t ${task.cpus} ref.fa ${inputReadFastq} > ${id}.paf
    """
}

process catReads {
    input:
    file('*.fastq.gz') from readsToCat.collect()

    output:
    file('all_reads.fastq.gz') into cattedReads

    """
    zcat *.fastq.gz | gzip > all_reads.fastq.gz
    """
}

process joinPafs {
    cpus 16
    publishDir 'aligned'

    input:
    file '*.paf' from aligned.collect()

    output:
    file 'all.paf' into alignedAll

    """
    cat *.paf > all.paf
    """
}

process splitByContig {
    input:
    file 'ref.fa' from reference
    file 'ref.fa.fai' from referenceFaidx

    output:
    file 'contig_*.fa' into contigs

    """
    cut -f1 ref.fa.fai | while read ctg; do
        samtools faidx ref.fa \$ctg > contig_\${ctg}.fa
    done
    """
}

contigs
    .flatMap()
    .map { tuple ((it.baseName =~ /contig_(.*)/)[0][1], it) }
    .combine(alignedAll)
    .combine(cattedReads)
    .set { contigsPafReads }

process racon {
    cpus 16

    input:
    set val(ctg),
        file('contig.fa'),
        file('all.paf'),
        file('all.fastq.gz') from contigsPafReads

    output:
    file("${ctg}_corrected.fa") into corrected

    """
    awk -v ctg=$ctg '\$6 == ctg' all.paf > contig.paf
    racon -t ${task.cpus} all.fastq.gz contig.paf contig.fa \
        > ${ctg}_corrected.fa
    """
}

process catContigs {
    publishDir 'out'

    input:
    file "*_corrected.fa" from corrected.collect()

    output:
    file "corrected.fa" into correctedCatted

    """
    cat *_corrected.fa > corrected.fa
    """
}


