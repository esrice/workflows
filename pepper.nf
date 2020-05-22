#!/usr/bin/env nextflow

params.reference = '../contigs/jackson-trio/jackson-haplotypepaternal/jackson-haplotypepaternal.contigs.fasta'
params.reads = '../contigs/jackson-trio/haplotype/haplotype-paternal.fasta.gz'

reference = file(params.reference)
Channel.fromPath(params.reads)
    .map { f -> tuple(f.baseName, f) }
    .set { readsWithId }

process align {
    cpus 28
    memory '230 GB'

    input:
    file 'ref.fa' from reference
    set id, file(inputReadFastq) from readsWithId

    output:
    file "${id}.bam" into aligned

    """
    minimap2 -ax map-ont -t ${task.cpus} ref.fa ${inputReadFastq} |
        samtools view -hb -q 60 -F 0x904 - |
        samtools sort -@ ${task.cpus} - > ${id}.bam
    """
}

process mergeBams {
    memory '128 GB'
    publishDir 'alignments'

    input:
    file bams from aligned.collect()

    output:
    file "merged.bam" into merged
    file "merged.bam.bai" into mergedIndex

    """
    bams="${bams}"
    words=( \$bams )
    if [[ "1" == "\${#words[@]}" ]]; then
        cp ${bams} merged.bam
    else
        samtools merge merged.bam ${bams}
    fi
    samtools index merged.bam
    """
}

process pepper {
    cpus 56
    memory '0'
    container 'kishwars/pepper:latest'
    publishDir 'pepper_out'
    clusterOptions = '--exclusive --account=johnsonlab'

    input:
    file "contigs.fa" from reference
    file "merged.bam" from merged
    file "merged.bam.bai" from mergedIndex

    output:
    file('pepper_out') into pepperOut

    """
    mkdir pepper_out
    pepper download_models -o models
    pepper polish -b merged.bam -f contigs.fa \
        -m models/pepper_r941_guppy344_human.pkl \
        -o pepper_out/pepper_out -tpc 7 -c 8 -t ${task.cpus}
    """
}

