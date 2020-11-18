from snakemake.utils import R
from os.path import join
import re,os

from os import listdir


#configfile: "CutAndRunConfig.json"

configfile: "run.json"

workpath = config['project']['workpath']
exp_name= config['project']['experiment_name']
pfamily=config['project']



# create output directories

fastq_dir='FASTQ'
trim_dir='trim'
bam_dir='bam'
split_bam_dir='split_bam'
bam_dir2='bam.120bp'
bigwig_dir='bigwig'
deeptools_dir='deeptools'
macs_dir='macs_peaks'
macs_dir2='macs_peaks.120bp'
seacr_dir='seacr_peaks'
seacr_dir2='seacr_peaks.120bp'
frip_dir='FRiP'

for d in [trim_dir,bam_dir,bam_dir2,split_bam_dir,bigwig_dir,deeptools_dir,macs_dir,seacr_dir,macs_dir2,seacr_dir2,frip_dir]:
	if not os.path.exists(join(workpath,d)):
		os.mkdir(join(workpath,d))


# define samples inputs and controls
samples, = glob_wildcards(join(fastq_dir, '{sample}_R1.fastq.gz'))

chips = config['project']['chips']
inputs = config['project']['controls']
chip2input = {} 
for chip in chips: 
    for input in inputs: 
        chip2input[chip] = input 



# rules

rule all:
     input:
        expand(join(workpath,"QC","{name}.qcmetrics"), name=samples),
        expand(join(workpath,"QC","{ext}.QCTable.txt"), ext=exp_name),
        expand(join(workpath,"multiqc_report.html")),
        expand(join(workpath,bigwig_dir,"{name}.mm10.RPKM.bw"),name=samples),
        expand(join(workpath,bigwig_dir,"{name}.mm10.120bp.RPKM.bw"),name=samples),
        expand(join(workpath,deeptools_dir,"spearman_heatmap.{ext}.pdf"),ext=exp_name),
        expand(join(workpath,deeptools_dir,"spearman_scatterplot.{ext}.pdf"),ext=exp_name),
        expand(join(workpath,deeptools_dir,"pca.{ext}.pdf"),ext=exp_name),
        expand(join(workpath,macs_dir,"{name}","{name}_peaks.narrowPeak"),name=chips),
        expand(join(workpath,seacr_dir,"{name}","{name}.seacr.norm.stringent.bed"),name=chips),
        expand(join(workpath,macs_dir2,"{name}","{name}_peaks.narrowPeak"),name=chips),
        expand(join(workpath,seacr_dir2,"{name}","{name}.seacr.norm.stringent.120bp.bed"),name=chips),
        expand(join(workpath,bam_dir2, "{name}.mm10.sorted.120bp.bam"), name=samples),
        expand(join(workpath,split_bam_dir,"{name}.mm10.sorted.pdf"), name=samples),
        expand(join(workpath,frip_dir,"{name}.FRiP_table.txt"),name=chips),
        expand(join(workpath,"QC","{ext}.FRiP_table.txt"),ext=exp_name),




rule trim_CutAndRun:
    input:
        file1=join(workpath,fastq_dir,"{name}_R1.fastq.gz"),
        file2=join(workpath,fastq_dir,"{name}_R2.fastq.gz"),
    output:
        outfq1=join(workpath,trim_dir,"{name}_R1.paired.fastq.gz"),
        outfq2=join(workpath,trim_dir,"{name}_R2.paired.fastq.gz"),
    params:
        rname="CR:trim",
        trimmomaticver=config['bin']['TRIMMOMATICVER'],
        workpath=config['project']['workpath'],
        fastawithadaptersetd=join(workpath,config['references']['FASTAWITHADAPTERSETD']),
        kseqbin=join(workpath,config['bin']['KSEQBIN']),
        minlen=25,
        leadingquality=20,
        trailingquality=20,
        javaram="64g",
        kseqlen=42,
    threads: 32
    shell: """
        module load {params.trimmomaticver};
        if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID ;fi
        cd /lscratch/$SLURM_JOBID;
        java -jar $TRIMMOJAR PE -threads {threads} -phred33 {input.file1} {input.file2}  R1.paired.fastq.gz R1.unpaired.fastq.gz R2.paired.fastq.gz R2.unpaired.fastq.gz ILLUMINACLIP:{params.fastawithadaptersetd}:2:15:4:4:true LEADING:{params.leadingquality} TRAILING:{params.trailingquality} SLIDINGWINDOW:4:15 MINLEN:{params.minlen};
        {params.kseqbin}/kseq_test R1.paired.fastq.gz {params.kseqlen} R1.step2.trim.fastq.gz;
        {params.kseqbin}/kseq_test R2.paired.fastq.gz {params.kseqlen} R2.step2.trim.fastq.gz;
        mv R1.step2.trim.fastq.gz {output.outfq1};
        mv R2.step2.trim.fastq.gz {output.outfq2};
        """

rule fastqc:
    input:
        expand(join(workpath,trim_dir,"{name}_R{rn}.paired.fastq.gz"), name=samples,rn=[1,2]),
    output: 
        join(workpath,"trim_QC"),
    params:
        fastqcver=config['bin']['FASTQCVER'],
        rname="fastqc",
    threads: 32
    shell:"""
    mkdir -p {output};
    module load {params.fastqcver};
    fastqc {input} -t {threads} -o {output}
    """

rule multiqc:
    input:
        join(workpath,"trim_QC"),
    params:
        rname="multiqc",
        multiqc=config['bin']['MULTIQCVER'],
    output:
        join(workpath,"multiqc_report.html")
    shell:"""
    module load {params.multiqc};
    multiqc {input}
    """



rule align_CutAndRun:
    input:
        file1=join(workpath,trim_dir,"{name}_R1.paired.fastq.gz"),
        file2=join(workpath,trim_dir,"{name}_R2.paired.fastq.gz"),
    params:
        d=join(workpath,bam_dir),
       	rname='CR:bowtie2',
        reference=join(workpath,config['references']['BOWTIE2']),
        bowtie2ver=config['bin']['BOWTIE2VER'],
        samtoolsver=config['bin']['SAMTOOLSVER'],
    output:
        outbam1=join(workpath,bam_dir,"{name}.sorted.bam"), 
        flagstat1=join(workpath,bam_dir,"{name}.sorted.bam.flagstat"),
    threads: 32
    shell: """ 
        module load {params.bowtie2ver};
        module load {params.samtoolsver};
        bowtie2 -p {threads} --dovetail --phred33  --very-sensitive -x {params.reference}  -1 {input.file1} -2 {input.file2}  | samtools view -bS - |  samtools sort -@{threads} -o {output.outbam1};
        samtools index {output.outbam1};
        samtools flagstat {output.outbam1} > {output.flagstat1};
        """

rule splitSpikein_CutAndRun:
    input:
        bam=join(workpath,bam_dir,"{name}.sorted.bam"),
    output:
        outbam1=join(workpath,split_bam_dir,"{name}.mm10.sorted.bam"),
        outbam2=join(workpath,split_bam_dir,"{name}.spikein.sorted.bam"),
        outbam3=join(workpath, bam_dir2,"{name}.mm10.sorted.120bp.bam"),
        outstats1=join(workpath, split_bam_dir,"{name}.mm10.sorted.bam.flagstat"),
        outstats2=join(workpath, split_bam_dir,"{name}.spikein.sorted.bam.flagstat"),
        outstats3=join(workpath, bam_dir2,"{name}.mm10.sorted.120bp.bam.flagstat"),
    params:
        d=join(workpath,bam_dir),
        rname='CR:split',
        chrs=config['references']['CHROMS'],
        spikeinchrs=config['references']['SPIKEINCHROMS'],
        samtoolsver=config['bin']['SAMTOOLSVER'],
        kseqbin=join(workpath,config['bin']['KSEQBIN']),
    threads: 32
    run:
        commoncmd="module load {params.samtoolsver}; "
        cmd1="samtools view -bh -@{threads} {input.bam} "+ " ".join(params.chrs)+" -o {output.outbam1}; samtools index {output.outbam1}; samtools flagstat {output.outbam1} > {output.outstats1}; "
        cmd2="samtools view -bh -@{threads} {input.bam} {params.spikeinchrs} -o {output.outbam2}; samtools index {output.outbam2}; samtools flagstat {output.outbam2} > {output.outstats2}; "
        cmd3="samtools view -h {output.outbam1} | LC_ALL=C awk -f {params.kseqbin}/filter_below.awk | samtools view -Sb - >{output.outbam3}; samtools index {output.outbam3}; samtools flagstat {output.outbam3} > {output.outstats3} "
        shell(commoncmd+cmd1+cmd2+cmd3)



rule bam2bigwig_CutAndRun:
     input:
        bam1=join(workpath,split_bam_dir,"{name}.mm10.sorted.bam"),
        bam2=join(workpath,bam_dir2,"{name}.mm10.sorted.120bp.bam"),
     output:
        bigwig1=join(workpath, bigwig_dir, "{name}.mm10.RPKM.bw"),
        bigwig2=join(workpath, bigwig_dir, "{name}.mm10.120bp.RPKM.bw"),
     params:
        rname="CR:bam2bw",
        deeptoolsver=config['bin']['DEEPTOOLSVER'],
     threads:32
     shell:"""
        module load {params.deeptoolsver};
        bamCoverage --bam {input.bam1} -o {output.bigwig1} --binSize 25 --smoothLength 75 --numberOfProcessors {threads} --normalizeUsing RPKM --centerReads --extendReads --ignoreForNormalization chrM;
        bamCoverage --bam {input.bam2} -o {output.bigwig2} --binSize 25 --smoothLength 75 --numberOfProcessors {threads} --normalizeUsing RPKM --centerReads --extendReads --ignoreForNormalization chrM;
        """

rule ppqt:
    input:
        bam=join(workpath, split_bam_dir,"{name}.mm10.sorted.bam"),
    output:
        ppqt=join(workpath, split_bam_dir,"{name}.mm10.sorted.ppqt"),
        pdf=join(workpath, split_bam_dir,"{name}.mm10.sorted.pdf"),
    params:
        rname="CR:ppqt",
        samtoolsver=config['bin']['SAMTOOLSVER'],
        rver=config['bin']['RVER'],
    run:
        commoncmd="module load {params.samtoolsver};module load {params.rver};"
        cmd="samtools view -b -f 66 -o /lscratch/$SLURM_JOBID/bam1.f66.bam {input.bam}; \
                samtools index /lscratch/$SLURM_JOBID/bam1.f66.bam; \
                Rscript Scripts/phantompeakqualtools/run_spp.R \
                -c=/lscratch/$SLURM_JOBID/bam1.f66.bam -savp={output.pdf} -out={output.ppqt} \
                -tmpdir=/lscratch/$SLURM_JOBID -rf;"
        shell(commoncmd+cmd)

rule NRF:
    input:
        bam=join(workpath,bam_dir,"{name}.sorted.bam"),
    params:
        rname='pl:NRF',
        samtoolsver=config['bin']['SAMTOOLSVER'],
        rver=config['bin']['RVER'],
        preseqver=config['bin']['PRESEQVER'],
        nrfscript=join(workpath,"Scripts","atac_nrf.py"),         
    output:
        preseq=join(workpath,"QC","{name}.preseq.dat"),
        preseqlog=join(workpath,"QC","{name}.preseq.log"),
        nrf=join(workpath,"QC","{name}.nrf"),
    threads: 16
    shell: """
module load {params.preseqver};
preseq lc_extrap -P -B -D -o {output.preseq} {input.bam} -seed 12345 -v -l 100000000000 2> {output.preseqlog}
python {params.nrfscript} {output.preseqlog} > {output.nrf}
        """


rule readstats:
    input:
        flagstat=join(workpath,bam_dir,"{name}.sorted.bam.flagstat"),
        infq=join(workpath,fastq_dir,"{name}_R1.fastq.gz"),
        infq2=join(workpath,trim_dir,"{name}_R1.paired.fastq.gz"),
        refflagstat=join(workpath,split_bam_dir,"{name}.mm10.sorted.bam.flagstat"),
        refflagstat2=join(workpath,bam_dir2,"{name}.mm10.sorted.120bp.bam.flagstat"),
        spikeinflagstat=join(workpath,split_bam_dir,"{name}.spikein.sorted.bam.flagstat"),
        nrf=join(workpath,"QC","{name}.nrf"),
    params:
        rname='CR:QCstats',
        filterCollate='Scripts/filterMetrics',   
    output:
        sampleQCfile=join(workpath,"QC","{name}.qcmetrics"),
    threads: 16
    shell: """
# Number of reads
zcat {input.infq} | wc -l | {params.filterCollate} {wildcards.name} tnreads > {output.sampleQCfile}
zcat {input.infq2} | wc -l | {params.filterCollate} {wildcards.name} tnreads >> {output.sampleQCfile}
# Number of mapped reads
grep 'mapped (' {input.flagstat} | awk '{{print $1,$3}}' | {params.filterCollate} {wildcards.name} mnreads >> {output.sampleQCfile}
# Number of ref reads
grep 'mapped (' {input.refflagstat} | awk '{{print $1,$3}}' | {params.filterCollate} {wildcards.name} mnreads >> {output.sampleQCfile}
grep 'mapped (' {input.refflagstat2} | awk '{{print $1,$3}}' | {params.filterCollate} {wildcards.name} mnreads >> {output.sampleQCfile}
# Number of spikein reads
grep 'mapped (' {input.spikeinflagstat} | awk '{{print $1,$3}}' | {params.filterCollate} {wildcards.name} mnreads >> {output.sampleQCfile}
# NRF, PBC1, PBC2
cat {input.nrf} | {params.filterCollate} {wildcards.name} nrf >> {output.sampleQCfile}
"""



rule readstats_table:
    input:
        expand(join(workpath,"QC","{name}.qcmetrics"), name=samples),
    params:
        rname='CR:QCTable',
        inputstring=",".join(expand(join(workpath,"QC","{name}.qcmetrics"), name=samples)),
        summary_script=join(workpath,"Scripts/compile_cut_and_run_stats.v2.pl")
    output:
        qctable=join(workpath,"QC","{exp_name}.QCTable.txt"),
    threads: 16
    shell: """
        perl {params.summary_script} {params.inputstring} > {output.qctable}
        """


rule QC_plots:
    input:
        expand(join(workpath,bigwig_dir,"{name}.mm10.RPKM.bw"),name=samples),
    output:
        heatmap=join(workpath,deeptools_dir,"spearman_heatmap.{exp_name}.pdf"),
        scatter=join(workpath,deeptools_dir,"spearman_scatterplot.{exp_name}.pdf"),
        pca=join(workpath,deeptools_dir,"pca.{exp_name}.pdf"),
        npz=temp(join(workpath,deeptools_dir,"{exp_name}.npz")),
    threads: 32
    params:
        rname="CR:deeptools_QC",
        deeptoolsver=config['bin']['DEEPTOOLSVER'],
    shell:"""
        module load {params.deeptoolsver};
        multiBigwigSummary bins -p {threads} -b {input}  --smartLabels -out {output.npz}
        plotCorrelation -in {output.npz} -o {output.heatmap} -c 'spearman' -p 'heatmap' --skipZeros --removeOutliers --plotNumbers
        plotCorrelation -in {output.npz} -o {output.scatter} -c 'spearman' -p 'scatterplot' --skipZeros --removeOutliers
        plotPCA -in {output.npz} -o {output.pca}
        """

### peak caling

rule macs_peakcall:
    input:
        chip = join(workpath,split_bam_dir,"{name}.mm10.sorted.bam"),
    output:
        join(workpath,macs_dir,"{name}","{name}_peaks.narrowPeak"),
    params:
        rname='CR:MACS',
        macsver=config['bin']['MACS2VER'],
    shell: """
        module load {params.macsver};
        macs2 callpeak -t {input.chip} -g mm -n {wildcards.name} --outdir {workpath}/{macs_dir}/{wildcards.name} -q 0.01 --keep-dup="all" -f "BAMPE";
    """


rule macs_peakcall_120:
    input:
        chip2 = join(workpath,bam_dir2,"{name}.mm10.sorted.120bp.bam"),
    output:
        join(workpath,macs_dir2,"{name}","{name}_peaks.narrowPeak"),
    params:
        rname='CR:MACS120',
        macsver=config['bin']['MACS2VER'],
    shell: """
        module load {params.macsver};
        macs2 callpeak -t {input.chip2} -g mm -n {wildcards.name} --outdir {workpath}/{macs_dir2}/{wildcards.name} -q 0.01 --keep-dup="all" -f "BAMPE";
    """


rule seacr_peakcall:
    input:
        chip = join(workpath,split_bam_dir,"{name}.mm10.sorted.bam"),
    output:
        stringent=join(workpath,seacr_dir,"{name}","{name}.seacr.norm.stringent.bed"),
        stringent_sum=join(workpath,seacr_dir,"{name}","{name}.seacr.norm.stringent.summits.bed"),
        relaxed=join(workpath,seacr_dir,"{name}","{name}.seacr.norm.relaxed.bed"),
        relaxed_sum=join(workpath,seacr_dir,"{name}","{name}.seacr.norm.relaxed.summits.bed"),
    params:
        rname='CR:seacr',
        macsver=config['bin']['MACS2VER'],
        bedtoolsver=config['bin']['BEDTOOLSVER'],
        samtoolsver=config['bin']['SAMTOOLSVER'],
        kseqbin=join(workpath,config['bin']['KSEQBIN']),
        seacrver=config['bin']['SEACRVER'],
        bedopsver=config['bin']['BEDOPSVER'],
        ctrl = lambda w : join(workpath,split_bam_dir,chip2input[w.name] + ".mm10.sorted.bam"),
    shell:"""
        module load {params.macsver}
        module load {params.bedtoolsver}
        module load {params.samtoolsver}
        if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID; fi 

        cd /lscratch/$SLURM_JOBID

        macs2 callpeak -t {input.chip} -g mm -f BAMPE -n treat --outdir tmp -q 0.01 -B --keep-dup all
        macs2 callpeak -t {params.ctrl} -g mm -f BAMPE -n ctrl --outdir tmp2 -q 0.01 -B --keep-dup all

        module load python/2.7
        python {params.kseqbin}/change.bdg.py tmp/treat_treat_pileup.bdg >treat.bdg
        python {params.kseqbin}/change.bdg.py tmp2/ctrl_treat_pileup.bdg >ctrl.bdg

        module load {params.seacrver}
        SEACR.sh treat.bdg ctrl.bdg norm stringent st.out
        SEACR.sh treat.bdg ctrl.bdg norm relaxed rel.out

        module load {params.bedopsver}
        sort-bed st.out.stringent.bed  >{output.stringent}
        sort-bed rel.out.relaxed.bed  >{output.relaxed}

        python {params.kseqbin}/get_summits_seacr.py {output.stringent} | sort-bed - >{output.stringent_sum}
        python {params.kseqbin}/get_summits_seacr.py {output.relaxed} | sort-bed -  >{output.relaxed_sum}
    """



rule seacr_peakcall_120:
    input:
        chip = join(workpath,bam_dir2,"{name}.mm10.sorted.120bp.bam"),
    output:
        stringent=join(workpath,seacr_dir2,"{name}","{name}.seacr.norm.stringent.120bp.bed"),
        stringent_sum=join(workpath,seacr_dir2,"{name}","{name}.seacr.norm.stringent.summits.120bp.bed"),
        relaxed=join(workpath,seacr_dir2,"{name}","{name}.seacr.norm.relaxed.120bp.bed"),
        relaxed_sum=join(workpath,seacr_dir2,"{name}","{name}.seacr.norm.relaxed.summits.120bp.bed"),
    params:
        rname='CR:seacr120',
        macsver=config['bin']['MACS2VER'],
        bedtoolsver=config['bin']['BEDTOOLSVER'],
        samtoolsver=config['bin']['SAMTOOLSVER'],
        kseqbin=join(workpath,config['bin']['KSEQBIN']),
        seacrver=config['bin']['SEACRVER'],
        bedopsver=config['bin']['BEDOPSVER'],
        ctrl = lambda w : join(workpath,bam_dir2,chip2input[w.name] + ".mm10.sorted.120bp.bam"),
    shell:"""
        module load {params.macsver}
        module load {params.bedtoolsver}
        module load {params.samtoolsver}
        if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID; fi 

        cd /lscratch/$SLURM_JOBID

        macs2 callpeak -t {input.chip} -g mm -f BAMPE -n treat --outdir tmp -q 0.01 -B --keep-dup all
        macs2 callpeak -t {params.ctrl} -g mm -f BAMPE -n ctrl --outdir tmp2 -q 0.01 -B --keep-dup all

        module load python/2.7
        python {params.kseqbin}/change.bdg.py tmp/treat_treat_pileup.bdg >treat.bdg
        python {params.kseqbin}/change.bdg.py tmp2/ctrl_treat_pileup.bdg >ctrl.bdg

        module load {params.seacrver}
        SEACR.sh treat.bdg ctrl.bdg norm stringent st.out
        SEACR.sh treat.bdg ctrl.bdg norm relaxed rel.out

        module load {params.bedopsver}
        sort-bed st.out.stringent.bed  >{output.stringent}
        sort-bed rel.out.relaxed.bed  >{output.relaxed}

        python {params.kseqbin}/get_summits_seacr.py {output.stringent} | sort-bed - >{output.stringent_sum}
        python {params.kseqbin}/get_summits_seacr.py {output.relaxed} | sort-bed -  >{output.relaxed_sum}
    """



 #       expand(join(workpath,frip_dir,"{name}_FRiP_barplot.png"),name=chips)

rule FRiP:
     input:
        bam =join(workpath,split_bam_dir,"{name}.mm10.sorted.bam"),
        bed1=join(workpath,seacr_dir,"{name}","{name}.seacr.norm.relaxed.bed"),
        bed2=join(workpath,macs_dir,"{name}","{name}_peaks.narrowPeak"),
     output:
        join(workpath,frip_dir,"{name}.FRiP_table.txt"),
     params:
        rname="CR:frip",
        pythonver="python/3.5",
        script=join(workpath,"Scripts","compute-frip.py"),
        genome = join(workpath,config['references']['REFLEN']),
        outroot=join(workpath,frip_dir,"{name}"),
     shell: """
module load {params.pythonver}
python {params.script} -p "{input.bed1} {input.bed2}" -b "{input.bam}" -g {params.genome} -o {params.outroot}
"""


rule compile_frip:
    input:
        expand(join(workpath,frip_dir,"{name}.FRiP_table.txt"), name=chips),
    params:
        rname='CR:CompileFRiP',
        inputstring=" ".join(expand(join(workpath,frip_dir,"{name}.FRiP_table.txt"), name=chips)),
    output:
        qctable=join(workpath,"QC","{exp_name}.FRiP_table.txt"),
    shell: """
        echo "bedtool,bedsample,bamsample,bamcondition,n_reads,n_overlap_reads,FRiP,n_basesM" | tr ',' '\t' > {output.qctable}
        cat {params.inputstring} | grep -v bedtool >> {output.qctable}
        """

