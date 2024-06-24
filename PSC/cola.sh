#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -t 1-200
#$ -notify
#$ -S /bin/bash
#$ -pe smp 8
##$ -M gianpiero.canessa@uai.cl
##$ -m e
#$ -N psc_1000

python bab2.py $SGE_TASK_ID 1000 1 50 
#python sims_fijo.py
