process FASTQC {
    clusterOptions "--nodes=3"
    publishDir "$params.outDir/FastQC/", mode:'copy'
    conda '/home/ndo/nextflow/envs/ngsArt.yml'

    input:
    tuple val (read_id), path(reads) 
    
    output:
    path("fastqc_${read_id}")
    
    script:
    """
    mkdir fastqc_${read_id}
    fastqc ${reads} -o fastqc_${read_id} -t 2
    """
}


process FASTP { 
    label 'short'
    publishDir "$params.outDir/Fastp/", mode:'copy'
    conda '/home/ndo/nextflow/envs/ngsArt.yml'

    input:
    tuple val(read_id), path(reads)
    
    output:
    tuple val(read_id), path("${read_id}_trimmed.fq.gz"), emit: trimmed
    path("${read_id}.fastp.json"), emit: json
    
    script:
    """
    fastp -i ${reads} -o ${read_id}_trimmed.fq.gz \
    -q ${params.minQual} -l ${params.minLen} -y --thread 2 --adapter_sequence=CTGTCTCTTATACACATCT --json ${read_id}.fastp.json
    """
}

process FASTQC_POST_TRIMMING {
    clusterOptions "--nodes=3"
    publishDir "$params.outDir/FastQC_Post_Trimming/", mode:'copy'
    conda '/home/ndo/nextflow/envs/ngsArt.yml'

    input:
    tuple val (read_id), path(reads) 
    
    output:
    path("fastqc_${read_id}")
    
    script:
    """
    mkdir fastqc_${read_id}
    fastqc ${reads} -o fastqc_${read_id} -t 2
    """
}
// Indexing reference genome
process REF_INDEX {
    label 'short'
    conda '/home/ndo/nextflow/envs/ngsArt.yml'

    input: 
    path 'ref'

    script:
    """
    bwa index ${ref}  
    """
}

process MAPPING {
    clusterOptions "--nodes=3"
    label 'high'
    publishDir "$params.outDir/BWA/", mode:'copy'
    conda '/home/ndo/nextflow/envs/ngsArt.yml'
    
    input:
    tuple val(read_id), path(reads)
    path ("ref/*")

    output:
    tuple val(read_id), path("${read_id}_raw.bam")

    script:
    """
    bwa mem -t 8 ref/*fasta ${reads} \
    | samtools view -uh -o ${read_id}_raw.bam
    """
}

process QUERY_SORT {
    label 'high'
    publishDir "$params.outDir/QuerySort/", mode:'copy'

    input:
    tuple val(read_id), path(reads)

    output:
    tuple val(read_id), path("${read_id}_query_sorted.bam")

    script:
    """
    gatk SortSam \
        --INPUT ${reads}\
        --OUTPUT ${read_id}_query_sorted.bam \
        --TMP_DIR ${params.tmp_dir} \
        --SORT_ORDER queryname
    """
}

process MARK_REMOVE_DUP {
    label 'norm'
    publishDir "$params.outDir/RemoveDup/", mode:'copy'
  

    input:
    tuple val(read_id), path(reads)

    output:
    tuple val(read_id), path("${read_id}_remove_dup.bam")

    script:
    """
    gatk MarkDuplicates\
      -I ${reads} \
      -O ${read_id}_remove_dup.bam \
      -M ${params.remove_dup_metrics} \
      -REMOVE_DUPLICATES true
    """
}

process COOR_SORT{
    label 'high'
    publishDir "$params.outDir/CoorSort/", mode:'copy'


    input:
    tuple val(read_id), path(reads)

    output:
    tuple val(read_id), path("${read_id}_coord_sort.bam")

    script:
    """
    gatk --java-options "-Xmx24G" SortSam \
        --INPUT ${reads}\
        --OUTPUT ${read_id}_coord_sort.bam \
        --TMP_DIR ${params.tmp_dir} \
        --SORT_ORDER coordinate
    """
}

process INTERSECT {
    label 'norm'
    publishDir "$params.outDir/Intersect/", mode:'copy'
    conda '/home/ndo/nextflow/envs/ngsArt.yml'
    
    input: 
    tuple val(read_id), path(reads)
    path giab_bed

    output:
    tuple val(read_id), path("${read_id}_intersected.bam")

    script:
    """
    bedtools intersect -header -a ${reads} -b ${giab_bed} > ${read_id}_intersected.bam
    samtools index ${read_id}_intersected.bam ${read_id}_intersected.bam.bai
    """
}

// Indexing final BAM file

process INDEXING_BAM {
    label 'norm'
    publishDir "$params.outDir/Intersect/", mode:'copy'
    conda '/home/ndo/nextflow/envs/ngsArt.yml'
    
    input: 
    tuple val(read_id), path(reads)

    output:
    tuple val(read_id), path("${read_id}_intersected.bam.bai")

    script:
    """
    samtools index ${read_id}_intersected.bam ${read_id}_intersected.bam.bai
    """
}


process VAR_CALLING {
    label 'high'
    clusterOptions "--nodes=2"
    publishDir "$params.outDir/VarCall/", mode:'copy'
    conda '/home/ndo/nextflow/envs/ngsArt.yml'

    input:
    tuple val(read_id), path(reads)
    path ("ref/*")

    output:
    tuple val(read_id), path("${read_id}_vardict.vcf")

    script:
    """
    vardict-java \
    -th 32 \
    -G ref/*fasta \
    -f ${params.AF_THR} \
    -b ${reads} \
    -c 1 -S 2 -E 3 -g 4 ${params.giab_bed} \
    | teststrandbias.R | var2vcf_valid.pl \
    -E -f ${params.AF_THR} > ${read_id}_vardict.vcf
    """
}

process BAM_TO_BED {
    label 'norm'
    clusterOptions "--nodes=2"
    publishDir "$params.outDir/SampleBed/", mode:'copy'
    conda '/home/ndo/nextflow/envs/ngsArt.yml'

    input:
    tuple val(read_id), path(reads)
   
    output:
    tuple val(read_id), path("${read_id}_bam.bed")

    script:
    """
    bedtools bamtobed -i ${reads} \
    > ${read_id}_bam.bed
    """
}

workflow {
  reads_ch = Channel
          .fromPath(params.raw_fastq_SE, checkIfExists: true)
          .map { tuple( it.baseName, it ) }
  FASTQC(reads_ch)
  FASTP(reads_ch)
  FASTQC_POST_TRIMMING(FASTP.out.trimmed)
  ref_ch = Channel.fromPath(params.refGenome)
  ref_files = file("${params.ref_dir}/*")
  MAPPING(FASTP.out.trimmed, ref_files)
  QUERY_SORT(MAPPING.out)
  MARK_REMOVE_DUP(QUERY_SORT.out)
  COOR_SORT(MARK_REMOVE_DUP.out)
  INTERSECT(COOR_SORT.out, params.giab_bed)
  INDEXING_BAM(INTERSECT.out)
  VAR_CALLING(INTERSECT.out, ref_files)
  BAM_TO_BED(INTERSECT.out)
}
