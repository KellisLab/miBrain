#!/usr/bin/env bash
RNA_PATH=$(readlink -f "$1")
NPROC=$(nproc)
GENOME_DIR="/net/bmc-lab5/data/kellis/group/Benjamin/ref/STAR_gencode43/"
GTF="/net/bmc-lab5/data/kellis/group/Benjamin/ref/gencode.v43.annotation.gtf"
type -P STAR >/dev/null || module load star/2.7.9a
type -P featureCounts >/dev/null || module load subread/1.6.2
find "${RNA_PATH}" -mindepth 2 -maxdepth 2 -type f -name '*.fastq' -exec dirname {} \; | sort | uniq | xargs -n1 -I{} sh -c "STAR --runThreadN ${NPROC} --runMode alignReads --quantMode GeneCounts --outSAMtype BAM Unsorted SortedByCoordinate --readFilesIn {}/*.fastq --outFileNamePrefix {}/ --genomeDir ${GENOME_DIR}"
find "${RNA_PATH}" -mindepth 2 -maxdepth 2 -type f -name '*.fastq.gz' -exec dirname {} \; | sort | uniq | xargs -n1 -I{} sh -c "STAR --runThreadN ${NPROC} --runMode alignReads --readFilesCommand zcat --quantMode GeneCounts --outSAMtype BAM Unsorted SortedByCoordinate --readFilesIn {}/*.fastq.gz --outFileNamePrefix {}/ --genomeDir ${GENOME_DIR}"
find "${RNA_PATH}" -mindepth 2 -maxdepth 2 -type f -name "Aligned.out.bam" -exec dirname {} \; | xargs -P0 -I{} sh -c "featureCounts -a ${GTF} -o {}/featureCounts {}/Aligned.out.bam"
