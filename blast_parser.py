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

def collect_first_records(blast_frame):
    blast_file = open(blast_frame,'r')
    Lines = blast_file.readlines()

    hash_start=False
    acc_dict = {}

    for line in Lines:
        if line[0]=="#":
            hash_start=True
            continue
        elif hash_start==True and line[0]!="#":
            acc = str.split(line, ",")[1]
            if acc not in acc_dict.keys():
                acc_dict[acc] = 1
            else:
                acc_dict[acc] = acc_dict[acc] + 1
        elif hash_start == True and line[0]=="#":
            hash_start=False
            continue

    acc_dict_sorted = sorted(acc_dict.items(), key=lambda x: x[1],reverse=True)       
    return acc_dict_sorted

def initiate_database():
    database = pd.DataFrame(['Accession','TaxID','Scientific'])
    return database

def query_ncbi(accesion_list,empty_datab):
    sci_dict = {}
    count = 0
    for item in accesion_list:
        acc = item[0]
        tax_handle = Entrez.esummary(db="Nucleotide", id=acc, retmode="xml")
        tax_records = Entrez.parse(tax_handle)
        for tax in tax_records:
            tax_id = int(tax["TaxId"])
            time.sleep(2/100)
            name_handle = Entrez.esummary(db="Taxonomy", id=tax_id, retmode="xml")
            name_records = Entrez.parse(name_handle)
            for name in name_records:
                curr_name = name['ScientificName']
                if curr_name not in sci_dict:
                    sci_dict[curr_name] = item[1]
                else:
                    sci_dict[curr_name] = sci_dict[curr_name] + item[1]
                time.sleep(2/100)
            name_handle.close()
        tax_handle.close()
        count = count + 1 
        print(count)
        

        #TODO rebuild using tax_handle = Entrez.efetch(db="Nucleotide", id=acc, retmode="xml")
        # record = Entrez.read(tax_handle)
        # record[0]['GBSeq_source']
        
    return sci_dict


def main(arg):
    Entrez.email = os.environ.get("EMAIL_KEY")
    Entrez.api_key = os.environ.get("SECRET_KEY")

    empty_db = initiate_database()

    blast_name = arg.blast_output
    blast_dict = collect_first_records(blast_name)

    sci_name_dict = query_ncbi(blast_dict, empty_db)
    
if __name__ == "__main__":
    # Build Argument Parser in order to facilitate ease of use for user
    parser = argparse.ArgumentParser(
        description="Digests blast format 6 outputs and reports counts of each taxa identified")
    parser.add_argument('-b', '--blast_output', action='store', required=True,
                        help="name of the blast output file", dest='blast_output')
    parser.add_argument('-d', '--db', action='store', required=True,
                        help="name of the database from previous runs", dest='database')
    
    args = parser.parse_args()
    main(args)