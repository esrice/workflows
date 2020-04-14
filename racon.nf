#!/usr/bin/env nextflow

params.reference = '../contigs/jackson-trio/jackson-haplotypematernal/jackson-haplotypematernal.contigs.fasta'
params.reads = '../contigs/reads/long_reads/*.fastq.gz'
params.numChunks = 100

reference = file(params.reference)
Channel.fromPath(params.reads)
    .map { f -> tuple(f.baseName, f) }
    .set { reads }

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
    set id, file(inputReadFastq) from reads

    output:
    file "${id}.paf" into aligned

    """
    minimap2 -x ava-ont -t ${task.cpus} ref.fa ${inputReadFastq} > ${id}.paf
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
    file 'all.paf' from alignedAll
    file 'ref.fa' from reference
    file 'ref.fa.fai' from referenceFaidx

    output:
    file 'contig_*.fa' into contigs
    file 'contig_*.paf' into pafs

    """
    cut -f1 ref.fa.fai | while read ctg; do
        samtools faidx ref.fa \$ctg > contig_\${ctg}.fa
        awk -v ctg=\$ctg '$6 == ctg' all.paf > contig_\${ctg}.paf
    done
    """
}

contigs
    .flatMap()
    .map { tuple ((it.baseName ~= /contig_(.*)/)[0][1], it) }
    .set { contigsIds }
pafs
    .flatMap()
    .map { tuple ((it.baseName ~= /contig_(.*)/)[0][1], it) }
    .set { pafsIds }
contigsIds.join(pafsIds).set { contigsPafs }

process racon {
    input:
    set val(ctg), file('contig.fa'), file('alignments.paf') from contigsPafs


}


