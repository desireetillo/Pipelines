#!/bin/bash

module load python/3.5
module load snakemake/5.1.3

snakemake all --delete-all-output
