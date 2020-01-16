#!/usr/bin/env nextflow

params.reference = '../../contigs/dam.arrow2.fasta'
params.classified_bam = '../read_classification/classified/maternal_only.bam'
params.enzyme = 'GATC'

reference_file = file(params.reference)
classified_bam = file(params.classified_bam)

process fasta_index {
    input:
    file reference from reference_file

    output:
    file "${reference}.fai" into faidx

    "samtools faidx ${reference}"
}

process filter_chimeras {
    cpus 2

    input:
    file 'in.bam' from classified_bam

    output:
    file 'forward.bam' into forwardBam
    file 'reverse.bam' into reverseBam

    """
    samtools view -bh -f64 in.bam | filter_chimeras.py - > forward.bam &
    samtools view -bh -f128 in.bam | filter_chimeras.py - > reverse.bam &
    wait
    """
}

process combine {
    publishDir 'alignments'

    input:
    file 'r1.bam' from forwardBam
    file 'r2.bam' from reverseBam

    output:
    file 'combined.bam' into combinedbam

    """
    combine_ends.py r1.bam r2.bam | samtools fixmate -m - - \
        | samtools sort - | samtools markdup -r - combined.bam
    """
}

process bam2bed {
    publishDir 'alignments'

    input:
    file 'combined.bam' from combinedbam

    output:
    file 'combined.bed' into combinedbed

    """
    bedtools bamtobed -i combined.bam | sort -k 4 > combined.bed
    """
}

process salsa {
    publishDir 'salsa_out'

    input:
    file 'ref.fa' from reference_file
    file 'ref.fa.fai' from faidx
    file 'combined.bed' from combinedbed

    output:
    file 'salsa/*' into salsa_out

    """
    python \$SALSA_DIR/run_pipeline.py -a ref.fa -l ref.fa.fai \
        -b combined.bed -e ${params.enzyme} -o salsa -m yes -i 10 -p yes
    """
}
