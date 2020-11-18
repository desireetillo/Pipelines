from os.path import join
import re,os

from os import listdir



configfile: "run.json"

workpath = config['project']['workpath']
exp_name= config['project']['experiment_name']


# create output directories

fastq_dir='FASTQ'
trim_dir='trim'
bam_dir='bam'
processed_dir='processed'
fastqc_dir='trim_QC'


# define samples inputs 
samples, = glob_wildcards(join(fastq_dir, '{sample}_R1.fastq.gz'))

for d in [trim_dir,bam_dir,processed_dir,fastqc_dir]:
    if not os.path.exists(join(workpath,d)):
        os.mkdir(join(workpath,d))


rule all:
    input:
        expand(join(workpath,"multiqc_report.html")),
        expand(join(workpath, processed_dir,"{name}_pe.CpG_report.txt.gz"),name=samples),
        expand(join(workpath,bam_dir, "{name}_pe.bam"),name=samples),
        expand(join(workpath,trim_dir,"{name}_R{rn}_val_{rn}.fq.gz"), name=samples,rn=[1,2]),


rule trim_rrbs:
    input:
        file1=join(workpath,fastq_dir,"{name}_R1.fastq.gz"),
        file2=join(workpath,fastq_dir,"{name}_R2.fastq.gz"),
    output:
        outfq1=join(workpath,trim_dir,"{name}_R1_val_1.fq.gz"),
        outfq2=join(workpath,trim_dir,"{name}_R2_val_2.fq.gz"),
    params:
        rname="RRBS:trim",
        trimgalorever=config['bin']['TRIMGALOREVER'],
        min_qual=20,
    threads: 8
    shell:"""
        module load {params.trimgalorever};
        trim_galore --quality {params.min_qual} \
        --phred33 --output_dir {trim_dir} \
        --gzip --rrbs --illumina \
        --paired --cores {threads} {input.file1} {input.file2}
      """


rule fastqc:
    input:
        expand(join(workpath,trim_dir,"{name}_R{rn}_val_{rn}.fq.gz"), name=samples,rn=[1,2]),
    output:
        expand(join(workpath,fastqc_dir,"{name}_R{rn}_val_{rn}_fastqc.html"), name=samples,rn=[1,2]),
    params:
        fastqcver=config['bin']['FASTQCVER'],
        rname="fastqc",
    threads: 32
    shell:"""
    mkdir -p {fastqc_dir};
    module load {params.fastqcver};
    fastqc {input} -t {threads} -o {fastqc_dir}
    """


rule multiqc:
    input:
        expand(join(workpath, fastqc_dir,"{name}_R{rn}_val_{rn}_fastqc.html"), name=samples,rn=[1,2]),
        expand(join(workpath, bam_dir, "{name}_pe.bam"), name=samples),
        expand(join(workpath, processed_dir,"{name}_pe.CpG_report.txt.gz"),name=samples),
    params:
        rname="multiqc",
        multiqc=config['bin']['MULTIQCVER'],
    output:
        join(workpath,"multiqc_report.html")
    shell:"""
    module load {params.multiqc};
    multiqc -f {trim_dir} {bam_dir} {fastqc_dir} {processed_dir}
    """

rule bismark_align:
    input:
        file1=join(workpath,trim_dir,"{name}_R1_val_1.fq.gz"),
        file2=join(workpath,trim_dir,"{name}_R2_val_2.fq.gz"),
    output:
        join(workpath,bam_dir, "{name}_pe.bam"),
    params:
        rname='RRBS:align',
        reference=config['references']['BISMARK'],
        bismarkver=config['bin']['BISMARKVER'],
    threads: 32
    shell:"""
    module load {params.bismarkver}
    bismark --phred33-quals --bowtie2 -p {threads} \
    --genome  {params.reference} \
    --nucleotide_coverage --output_dir {bam_dir} \
    --basename {wildcards.name} \
    --fastq --temp_dir /lscratch/$SLURM_JOBID -1 {input.file1} -2 {input.file2}
    """

rule methylation_extract:
    input:
        join(workpath, bam_dir, "{name}_pe.bam"),
    output:
        join(workpath, processed_dir,"{name}_pe.CpG_report.txt.gz"),
    params:
        rname='RRBS:methyl_extract',
        reference=config['references']['BISMARK'],
        bismarkver=config['bin']['BISMARKVER'],
    shell:"""
    module load {params.bismarkver};
    bismark_methylation_extractor -p --bedGraph --buffer_size 80% --ignore_r2 2  --comprehensive --merge_non_CpG \
    --cytosine_report --genome_folder {params.reference} --gzip --output {processed_dir} {input}
    """


