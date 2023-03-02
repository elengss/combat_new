#!/bin/bash
#$ -wd /well/combat/users/nif917/ref/CD8/logs
#$ -P combat.prjc
#$ -N my-python-job
#$ -q short.qc
#$ -pe shmem 4
cd /well/combat/users/nif917/ref/CD8/
module load R-bundle-Bioconductor
Rscript prep_covid_tem_RData.R
