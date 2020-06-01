#!/usr/bin/env nextflow

params.damReference = '../pepper_maternal/pepper_out/pepper_out/pepper_out_pepper.fa'
params.sireReference = '../pepper_paternal/pepper_out/pepper_out/pepper_out_pepper.fa'
params.hicReads = 'reads/HiC_*_R{1,2}_001.fastq.gz'
params.enzyme = 'GATC'

damReference = file(params.damReference)
sireReference = file(params.sireReference)
hicReads = Channel.fromFilePairs(params.hicReads)

process fasta_index {
    input:
    file 'dam.fa' from damReference
    file 'sire.fa' from sireReference

    output:
    file "dam.fa.fai" into damFaidx
    file "sire.fa.fai" into sireFaidx

    """
    samtools faidx dam.fa
    samtools faidx sire.fa
    """
}

process bwaIndex {
    module 'bwa'
    publishDir 'bwa_idx'
    memory '325 GB'
    cpus 2

    input:
    file 'dam.fa' from damReference
    file 'sire.fa' from sireReference

    output:
    file 'dam.fa*' into damIndex
    file 'sire.fa*' into sireIndex

    """
    bwa index dam.fa &
    bwa index sire.fa &
    wait
    """
}

hicReads.into { hicReadsDam; hicReadsSire }

process alignDam {
    module 'bwa/bwa-0.7.17:samtools/samtools-1.9'
    publishDir 'unclassified'
    cpus 14
    memory '125 GB'

    input:
    file idx from damIndex
    set( val (id), file (fastqs) ) from hicReadsDam

    output:
    set (val(id), file("dam.${id}.unclassified.bam")) into damBamsUnclassified

    """
    bwa mem -t ${task.cpus} dam.fa ${fastqs} | samtools view -bh - \
        | samtools sort -n - > dam.${id}.unclassified.bam
    """
}

process alignSire {
    module 'bwa/bwa-0.7.17:samtools/samtools-1.9'
    publishDir 'unclassified'
    cpus 14
    memory '125 GB'

    input:
    file idx from sireIndex
    set( val (id), file (fastqs)) from hicReadsSire

    output:
    set (val(id), file("sire.${id}.unclassified.bam")) into sireBamsUnclassified

    """
    bwa mem -t ${task.cpus} sire.fa ${fastqs} | samtools view -bh - \
        | samtools sort -n - > sire.${id}.unclassified.bam
    """
}

damBamsUnclassified.join(sireBamsUnclassified).set { damSireUnclassified }

process classify {
    publishDir 'classified'

    input:
    set(val(id), file(damBam), file(sireBam)) from damSireUnclassified

    output:
    tuple(val(id), file("dam.${id}.classified.bam")) into damClassified
    tuple(val(id), file("sire.${id}.classified.bam")) into sireClassified

    """
    classify_by_alignment --hapA-in ${damBam} \
        --hapA-out dam.${id}.classified.bam \
        --hapB-in ${sireBam} \
        --hapB-out sire.${id}.classified.bam
    """
}

damClassified
    .map { tuple('dam', it[0], it[1]) }
    .mix(sireClassified
        .map { tuple ('sire', it[0], it[1]) }
    ).set { allClassified }

process filter_chimeras {
    cpus 2
    container 'esrice/hic-pipeline:latest'

    input:
    tuple(val(parent), val(id), file("${parent}.${id}.bam")) from allClassified

    output:
    tuple(val(parent),
          val(id),
          file("forward.${parent}.${id}.bam"),
          file("reverse.${parent}.${id}.bam")) into filteredBams

    """
    samtools view -bh -f64 ${parent}.${id}.bam | filter_chimeras.py - \
        > forward.${parent}.${id}.bam &
    samtools view -bh -f128 ${parent}.${id}.bam | filter_chimeras.py - \
        > reverse.${parent}.${id}.bam &
    wait
    """
}

process combine {
    container 'esrice/hic-pipeline:latest'
    publishDir 'alignments'

    input:
    tuple(val(parent), val(id),
          file("r1.bam"),
          file("r2.bam")) from filteredBams

    output:
    tuple(val(parent), file("combined.${parent}.${id}.bed")) into combinedBeds

    """
    combine_ends.py r1.bam r2.bam | samtools fixmate -m - - \
        | samtools sort - | samtools markdup -r - combined.bam
    bedtools bamtobed -i combined.bam > combined.${parent}.${id}.bed
    """
}

combinedBeds
    .branch {
        dam: it[0] == 'dam'
        sire: it[0] == 'sire'
    }.set { combinedBranchedBeds }
combinedBranchedBeds.dam.map { it[1] }.set { damBeds }
combinedBranchedBeds.sire.map { it[1] }.set { sireBeds }

process salsaDam {
    container 'esrice/hic-pipeline:latest'
    publishDir 'salsa_out'

    input:
    file 'ref.fa' from damReference
    file 'ref.fa.fai' from damFaidx
    file beds from damBeds.collect()

    output:
    file 'salsa/*' into salsa_dam

    """
    cat ${beds} | sort -k 4 > all.bed
    python \$SALSA_DIR/run_pipeline.py -a ref.fa -l ref.fa.fai \
        -b all.bed -e ${params.enzyme} -o salsa -m yes -i 10 -p yes
    """
}

process salsaSire {
    container 'esrice/hic-pipeline:latest'
    publishDir 'salsa_out'

    input:
    file 'ref.fa' from sireReference
    file 'ref.fa.fai' from sireFaidx
    file beds from sireBeds.collect()

    output:
    file 'salsa/*' into salsa_sire

    """
    cat ${beds} | sort -k 4 > all.bed
    python \$SALSA_DIR/run_pipeline.py -a ref.fa -l ref.fa.fai \
        -b all.bed -e ${params.enzyme} -o salsa -m yes -i 10 -p yes
    """
}
