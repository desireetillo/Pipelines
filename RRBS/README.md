# RRBS pipelines

Methylation calling for RRBS data

These workflows process reads from RRBS experiments.  Reads are trimmed using trimmomatic with RRBS settings.  Trimmed reads are then mapped and methylation called using the bismark suite.  These workflows also  generates signal tracks (bigwigs) of methylation % and read coverage over CpG sites.  

**Note**:  may need to modify the *--ignore N* parameter for the bismark\_methlylation\_extractor step, depending on what the m-bias plots look like.


Two workflows currently:

* RRBS-SE.Snakefile - for single-read data
* RRBS-PE.Snakefile - for paired-end data 
*
Dependencies, locations of annotations are in the file run.json

A set of utility scripts in the "Scripts/" directory

## Setting up a run on biowulf

Need bismark formatted index (see script index_genome.sh for details)

Data should be in a directory FASTQ/  
Paired-end data should end with _R?.fastq.gz

Construct a meta.tab file (tab-delimited file containing IDs and treatment vector).  This file is ultimately used for running methylKit for the summary and DMR calling steps (currently run by hand, no rule yet).


Then run:

```
python Scripts/make_config.py \
--prefix prefix \
--template Templates/template_RRBS.json \
--meta meta.tab \
--exp_type RRBS
```

This generates the "run.json" file, which contains all of the parameters for the run

Do a dry run: 

```
module load snakemake/5.1.3
snakemake -n -s RRBS-SE.Snakefile
```

Submit job to biowulf:

```
sbatch run_analysis.sh
```

