# Jackwood Blast Parser
Script to Parse Blast type outputs to verify genome assemblies

## Install

```bash
$ source ~/.bashrc 
$ conda env create -f environment.yml
$ conda activate jackwood_blast
``` 

## Usage

```bash
$ python3 blast_parser.py -h
usage: blast_parser.py [-h] -b BLAST_OUTPUT [-d DATABASE] [-o DATABASE_FILE] [-n RESULT_FILE] [-v]

Digests blast format 6 outputs and reports counts of each taxa identified

optional arguments:
  -h, --help            show this help message and exit
  -b BLAST_OUTPUT, --blast_output BLAST_OUTPUT
                        name of the blast output file
  -d DATABASE, --db DATABASE
                        name of the database from previous runs
  -o DATABASE_FILE, --dbo DATABASE_FILE
                        name of the database to save to disk
  -n RESULT_FILE, --res_o RESULT_FILE
                        name of the result to save to disk
  -v, --verb            verbose output to the terminal
```

## Example Command

```bash
$ python3 blast_parser.py -b Advent_18s_2_result_problem_file_8_22.out -d database_1.csv -o database_2.csv -n Advent_18s_2_results.csv
```

## Getting the environment setup To use Docker
```bash
COMPUTER$ docker run -v $PWD:/mnt -it lorentzb/jackwood_blast

(base) root@c90fac8ae3ce:/build_conda_env# cd /mnt
(base) root@c90fac8ae3ce:/mnt# conda activate jackwood_blast
(jackwood_blast) root@c90fac8ae3ce:/mnt# source env_vars.sh
(jackwood_blast) root@c90fac8ae3ce:/mnt# python3 blast_parser.py -b BLAST_RESULTS.out -d DATABASE_IN.csv -o DATABASE_OUT.csv -n LABELD_RESULTS.csv

OPTIONAL 

$ mkdir result_files
$ mkdir databases

(jackwood_blast) root@c90fac8ae3ce:/mnt# python3 blast_parser.py -b BLAST_RESULTS.out -d databases/DATABASE_IN.csv -o databases/DATABASE_OUT.csv -n result_files/LABELD_RESULTS.csv
```