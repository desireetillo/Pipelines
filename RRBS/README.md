# RRBS pipelines

Methylation calling for RRBS data

This workflow takes fastq files, trims the data using trimmomatic with RRBS settings, maps the trimmed reads to the genome and calls methylation using the bismark suite,
and generates signal tracks (bigwigs).  

**Note**:  may need to modify the --ignore N parameter for the bismark_methlylation_extractor, depending on what the m-bias plots look like.


Two workflows currently:

RRBS-SE.Snakefile  - for single end data

RRBS-PE.Snakefile - for paired-end data (is missing rule to generate bigwigs).

Dependencies, locations of annotations are in the file run.json

A set of utility scripts in the "Scripts/" directory

## Setting up a run on biowulf

Need bismark formatted index (see script index_genome.sh for details)

Data in directory FASTQ/

Construct a meta.tab file (tab-delimited file containing IDs and treatment vector), ultimately used for methylKit QC as well as calling of DMRs (currently run by hand)



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
`module load snakemake/5.1.3;  
snakemake -n -s RRBS-SE.Snakefile`

To submit to biowulf:
`sbatch run_analysis.sh`
