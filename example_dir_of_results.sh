#!/bin/bash

#SBATCH --partition=batch
#SBATCH --job-name=blast_parse
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --mem=32gb

#SBATCH --mail-user=""
#SBATCH --mail-type=BEGIN,END,FAIL


eval "$(conda shell.bash hook)"
conda activate jackwood_blast

dt=$(date '+%d-&m-%Y_%H.%M.%S');

source env_vars.sh

blast_files=`ls /home/bj44874/Blast_results`
for FILE in $blast_files
do

        python3 blast_parser.py -b ../Blast_results/${FILE} -d database.csv -o database.csv -n ${FILE::-4}_$dt.csv

done
