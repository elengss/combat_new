#!/bin/bash
#$ -wd /well/combat/users/nif917/ref/CD8/logs
#$ -P combat.prjc
#$ -N my-python-job
#$ -q short.qc
cd /well/combat/users/nif917/ref/CD8/scripts_flu_beta_temra
module load R-bundle-Bioconductor
Rscript beta_zzz.R
