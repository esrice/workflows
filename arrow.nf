#!/usr/bin/env nextflow

params.reference = '../polishing_nf/arrow_out/arrow.fa'
params.subreads_folder = '../../reads/bams/'
params.numChunks = 100
params.bin = '~/htc/software/t/smrtlink/smrtcmds/bin'

reference_fasta = file(params.reference)
subreads = Channel.fromFilePairs(params.subreads_folder + '/*.{bam,bam.pbi}')

process make_reference {
    input:
    file 'ref.fa' from reference_fasta

    output:
    file 'ref.xml' into reference_xml

    """
    $params.bin/dataset create --type ReferenceSet --generateIndices \
        ref.xml ref.fa
    """
}

reference_xml.into { reference_xml_for_arrow; reference_xml_for_align }

process align {
    cpus 16

    input:
    file 'ref.xml' from reference_xml_for_align
    set id, file(subreads_bam_and_pbi) from subreads

    output:
    file "${id}.aln.bam" into aligned

    """
    TMPDIR=\$(pwd) $params.bin/pbmm2 align --sort -j 16 --preset SUBREAD \
        ref.xml ${id}.bam ${id}.aln.bam
    """
}

process split_by_contig {
    input:
    file '*.aln.bam' from aligned.collect()

    output:
    file 'all_aligned.*.xml' into aligned_split

    """
    $params.bin/dataset create --generateIndices --type AlignmentSet \
        all_aligned.xml *.aln.bam
    $params.bin/dataset split --contig \
        --chunks $params.numChunks all_aligned.xml
    """
}

process arrow {
    cpus 16

    input:
    file 'ref.xml' from reference_xml_for_arrow
    file aligned_chunk from aligned_split.flatten()

    output:
    file "${aligned_chunk.baseName}.fa" into fasta_chunks
    file "${aligned_chunk.baseName}.fq" into fastq_chunks

    """
    $params.bin/gcpp $aligned_chunk -j 16 -q 20 -r ref.xml \
        -o ${aligned_chunk.baseName}.fa,${aligned_chunk.baseName}.fq
    """
}

process merge {
    publishDir 'arrow_out'

    input:
    file "*.fa" from fasta_chunks.collect()
    file "*.fq" from fastq_chunks.collect()

    output:
    file "arrow.fa" into corrected_ref_fasta
    file "arrow.fq" into corrected_ref_fastq

    """
    cat *.fa > arrow.fa
    cat *.fq > arrow.fq
    """
}
