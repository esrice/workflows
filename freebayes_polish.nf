params.assembly = '/storage/htc/warrenlab/users/esrbhb/killifish/assembly/chrom_assignment/chroms.fa'
params.maxCoverage = 500

Channel.fromSRA('SRX4394280').set { short_reads }
Channel.fromPath("${params.assembly}.fai")
    .splitCsv(header: ['name', 'length', 'a', 'b', 'c'], sep: "\t")
    .set { regions }

process align {
    module 'bwa/bwa-0.7.17'
    cpus 26
    publishDir 'bams'
    memory '225 GB'

    input:
    set sra_id, file(reads) from short_reads

    output:
    file "${sra_id}.bam" into bam

    """
    bwa mem -t ${task.cpus} ${params.assembly} ${reads} | \
        samtools view -bh - | samtools fixmate -m - - | \
        samtools sort -T . -@ ${task.cpus} -m 2G - | \
        samtools markdup - ${sra_id}.bam
    """
}

process merge {
    input:
    file bams from bam.collect()

    output:
    file("merged.bam") into mergedBam

    """
    bams="${bams}"
    words=( \$bams )
    if [[ "1" == "\${#words[@]}" ]]; then
        cp ${bams} merged.bam
    else
        samtools merge merged.bam ${bams}
    fi
    """
}

mergedBam.combine(regions).set { bam_and_regions }

process freebayes {
    module 'bcftools/bcftools-1.8'

    input:
    set file("merged.bam"), val(region) from bam_and_regions

    output:
    file "${region.name}.bcf" into region_bcfs

    """
    samtools index merged.bam
    freebayes --bam merged.bam --region ${region.name}:1-${region.length} \
        --skip-coverage 500 -f ${params.assembly} \
        | bcftools view --no-version -Ob -o ${region.name}.bcf
    """
}

process concat_and_consensus {
    module 'bcftools/bcftools-1.8'
    publishDir 'consensus'

    input:
    file "*.bcf" from region_bcfs.collect()

    output:
    file "polished.fa" into polished
    file "joined.bcf" into variants
    file "report.txt" into report

    """
    bcftools concat -n *.bcf | bcftools view -Ou -e'type="ref"' \
        | bcftools norm -Ob -f ${params.assembly} -o joined.bcf
    bcftools index joined.bcf
    bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' \
        -Hla -f ${params.assembly} joined.bcf > polished.fa

    snp=\$(bcftools view -i'QUAL>1 && (GT="AA" || GT="Aa")' \
        -Ha joined.bcf | grep -c 'TYPE=snp')
    ins=\$(bcftools view -i'QUAL>1 && (GT="AA" || GT="Aa")' \
        -Ha joined.bcf | grep -c 'TYPE=ins')
    del=\$(bcftools view -i'QUAL>1 && (GT="AA" || GT="Aa")' \
        -Ha joined.bcf | grep -c 'TYPE=del')
    echo "\${snp} SNPs, \${ins} insertions, and \${del} deletions corrected." \
        > report.txt
    """
}

