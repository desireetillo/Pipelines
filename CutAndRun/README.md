# Cut&Run pipeline

## Overview

A workflow to process Cut&Run data on the [NIH biowulf cluster](https://hpc.nih.gov). Most steps and some borrowed code from this [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1802-4), as well as [CCBR's Pipeliner](https://github.com/CCBR/Pipeliner).


## Dependencies


	snakemake/5.1.3
	bedops/2.4.30  
	bedtools/2.27.1 
	bowtie/2-2.2.6  
	deeptools/3.0.1  
	fastqc/0.11.5  
	macs/2.1.1.20160309  
	multiqc/1.4  
	preseq/2.0.3  
	samtools/1.6  
	SEACR/1.3  
	trimmomatic/0.36  
	R/3.5  
	python/2.7 
	python/3.6 


## Setup

There are a set set of utility scripts in the "Scripts/" directory.  If pulling from github, must run the following before running the pipeline:

`chmod +x Scripts/CutAndRunTools/kseq_test Scripts/perl_lib/*pl`


Dependencies and paths reference annotation are set in the file Templates/template_CutAndRunConfig.json and may be edited to suit your purposes.

---


## Setting up a run on biowulf

Put fastqs in directory FASTQ/

read1 and read2 must end with _R?.gz

Construct a pairs.tab file (tab-delimited file containing IDs of IP and control sample), used for peak-calling:

	#IP<tab>control
	sample1<tab>control1
	sample2<tab>control2

Then run:

`python Scripts/make_config.py --prefix CutAndRun_Test --template Templates/template_CutAndRunConfig.json --pairs pairs.tab`

This generates the "run.json" file, which contains all of the parameters for the run.



To do a dry run: 

`module load snakemake/5.1.3`


`snakemake -n` 

To submit the workflow on biowulf (NIH HPC):

`sbatch run_analysis.sh`

---

## Details

1. **Trimming (trim_CutAndRun)**

    Tools:

    - trimmomatic/0.36
    - kseq\_test (if grabbing from github, do chmod +x Scripts/CutAndRunTools/kseq_test Scripts/perl\_lib/*pl)

    Trim sequences using `trimmomatic`

    ```
    java -jar $TRIMMOJAR PE -threads {threads} -phred33 {input.file1} {input.file2} R1.paired.fastq.gz R1.unpaired.fastq.gz R2.paired.fastq.gz R2.unpaired.fastq.gz ILLUMINACLIP:{params.fastawithadaptersetd}:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25;
    ```

    Then use another tool, `kseq` to  trim up to 6-bp adapters from the 3′ end of each read that was not effectively processed by Trimmomatic.

    ```
    {params.kseqbin}/kseq_test R1.paired.fastq.gz {params.kseqlen} R1.step2.trim.fastq.gz;
    {params.kseqbin}/kseq_test R2.paired.fastq.gz {params.kseqlen} R2.step2.trim.fastq.gz;
    ```

2. **Align using bowtie2 (align_CutAndRun)**

    Tools: 

    - bowtie/2-2.2.6

    ```
    bowtie2 -p {threads} --dovetail --phred33 --very-sensitive -x {params.reference} -1 {input.file1} -2 {input.file2}
    ```

    reference is a joint mm10-ecoli index (in db/bowtie2_index)

3.  **Filter alignments (splitSpikein_CutAndRun)**

    Tools:

    - samtools/1.6
    - filter_below.awk

    subset bamfiles into mm10 and E.coli alignments, compute stats using `samtools flagstat`

    extract ≤ 120-bp fragments from mm10 bamfile 

    Code for extracting 120bp fragments (used for peak-calling)

    ```
    samtools view -h {output.outbam1} | LC_ALL=C awk -f {params.kseqbin}/filter_below.awk | samtools view -Sb - >{output.outbam3}; samtools index {output.outbam3} 
    ```

4. **Generate signal track files from mm10 alignments (bam2bigwig_CutAndRun)**

tools: 

- deeptools/3.0.1

```
bamCoverage --bam {input.bam} -o {output.bigwig} --binSize 25 --smoothLength 75 --numberOfProcessors {threads} --normalizeUsing RPKM --centerReads --extendReads --ignoreForNormalization chrM;
```

5.**Call peaks on mm10 alignments (all alignments and ≤ 120bp fraction)**

rule **macs_peakcall**

tools:  

- macs/2.1.1.20160309

macs2 parameters (note that we do not use the control bam file for macs2, since it expects an input genomic control, whereas the typical control for Cut&Run experiments is an IgG control):

```
macs2 callpeak -t {ip.bam} -g mm -n {[wildcards.name](http://wildcards.name/)} --outdir {workpath}/{macs_dir2}/{[wildcards.name](http://wildcards.name/)} -q 0.01 --keep-dup="all" -f "BAMPE";
```

rule **seacr_peakcall:**

tools: 

- macs/2.1.1.20160309
- seacr/1.3
- bedops (but maybe can just use a sort call here)

code for seacr (see Snakefile). Requires conversion to bedgraph, and adjustment of coverage to integers 

```
# create coverage bedgraph file from bam file
macs2 callpeak -t {ip.bam} -g mm -f BAMPE -n treat --outdir tmp -q 0.01 -B --keep-dup all
macs2 callpeak -t {ctrl.bam} -g mm -f BAMPE -n ctrl --outdir tmp2 -q 0.01 -B --keep-dup all

# needs python/2.7
python {params.kseqbin}/change.bdg.py tmp/treat_treat_pileup.bdg >treat.bdg
python {params.kseqbin}/change.bdg.py tmp2/ctrl_treat_pileup.bdg >ctrl.bdg
SEACR.sh treat.bdg ctrl.bdg norm stringent st.out
SEACR.sh treat.bdg ctrl.bdg norm relaxed rel.out

# can probably use sort -k1,1 -k2,2n instead of using sort-bed from bedops
sort-bed st.out.stringent.bed  >{output.stringent}
sort-bed rel.out.relaxed.bed  >{output.relaxed}

python {params.kseqbin}/get_summits_seacr.py {output.stringent} | sort-bed - >{output.stringent_sum}
python {params.kseqbin}/get_summits_seacr.py {output.relaxed} | sort-bed -  >{output.relaxed_sum}
```

**QC steps:**

1. Run fastqc on trimmed sequences(**rule fastqc**; **rule multiqc**)
2. compute library complexity (**rule NRF**)
3. Compute read and alignment statistics (samtools flagstat on all bams, done in alignment step)
4. Deeptools QC plots (PCA, scatterplot, heatmap) from bigwigs (**rule QC_plots**)
5. FRiP score on peak calls (**rule FRiP**)


## **Outputs**


All processed files will be found in the following directories: 


**Trimmed fastqs:**

	trim/ 
	trim_QC/

**Alignments**: 

	bam/ 
	split_bam/
	bam.120bp/

**Signal tracks:**

	bigwig/

**QC:**

	FRiP/
	deeptools/
	QC/ ← read/alignment statistics here
	multiqc_data/

**Peak calls:**
	
	seacr_peaks/
	macs_peaks/
	macs_peaks.120bp/
	seacr_peaks.120bp/
