#!/usr/bin/env python3

import Bio
from Bio import Entrez
import time
import os
import argparse
import os
from os import path
import pandas as pd 
import numpy as np

def load_blast_table(blast_table):
    return []

def collect_first_records(blast_frame):
    return []



def main(arg):
    Entrez.email = os.environ.get("EMAIL_KEY")
    Entrez.api_key = os.environ.get("SECRET_KEY")

    blast_name = arg.blast_output
    blast_dataframe = load_blast_table(blast_name)
    
    collect_first_records(blast_dataframe)
    

if __name__ == "__main__":
    # Build Argument Parser in order to facilitate ease of use for user
    parser = argparse.ArgumentParser(
        description="Digests blast format 6 outputs and reports counts of each taxa identified")
    parser.add_argument('-b', '--blast_output', action='store', required=True,
                        help="name of the blast output file", dest='metadata')
    parser.add_argument('-d', '--db', action='store', required=True,
                        help="name of the database from previous runs", dest='ioi')
    
    args = parser.parse_args()
    main(args)