#!/usr/bin/env nextflow

params.reference = '../contigs/jackson-trio/jackson-haplotypepaternal/jackson-haplotypepaternal.contigs.fasta'
params.reads = '../contigs/reads/long_reads/*.fastq.gz'

reference = file(params.reference)
Channel.fromPath(params.reads)
    .map { f -> tuple(f.baseName, f) }
    .set { readsWithId }

process align {
    cpus 16
    memory '128 GB'

    input:
    file 'ref.fa' from reference
    set id, file(inputReadFastq) from readsWithId

    output:
    file "${id}.bam" into aligned

    """
    minimap2 -ax map-ont -t ${task.cpus} ref.fa ${inputReadFastq} |
        samtools view -hb -q 60 -F 0x904 - |
        samtools sort -@ ${task.cpus} - > ${id}.bam
    samtools index ${id}.bam
    """
}

process mergeBams {
    memory '128 GB'
    publishDir 'alignments'

    input:
    file "*.bam" from aligned.collect()

    output:
    file "merged.bam" into merged
    file "merged.bam.bai" into mergedIndex

    """
    samtools merge merged.bam *.bam
    samtools index merged.bam
    """
}

process marginpolish {
    cpus 56
    memory '0'
    container 'kishwars/helen:latest'
    publishDir 'mp_out'
    clusterOptions = '--exclusive --account=johnsonlab'

    input:
    file "contigs.fa" from reference
    file "merged.bam" from merged
    file "merged.bam.bai" from mergedIndex

    output:
    file('mp_out') into mpOut

    """
    mkdir mp_out
    helen download_models -o models
    marginpolish merged.bam contigs.fa models/MP_r941_guppy344_human.json \
        -o mp_out/mp_out -t ${task.cpus} -f
    """
}

process helen {
    cpus 16
    memory '200 GB'
    container 'kishwars/helen:latest'
    publishDir 'helen_out'

    input:
    file "mp_out" from mpOut

    output:
    file('helen_out') into corrected

    """
    helen download_models -o models
    helen polish -i mp_out/ -m models/HELEN_r941_guppy344_human.pkl \
        -t ${task.cpus} -o helen_out
    """
}

