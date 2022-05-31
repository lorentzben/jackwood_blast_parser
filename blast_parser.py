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
import progressbar
#import requests
#import lxml
#from bs4 import BeautifulSoup


def collect_first_records(blast_frame):
    blast_file = open(blast_frame, 'r')
    Lines = blast_file.readlines()

    hash_start = False
    acc_dict = {}

    for line in Lines:
        if line[0] == "#":
            hash_start = True
            continue
        elif hash_start == True and line[0] != "#":
            acc = str.split(line, ",")[1]
            if acc not in acc_dict.keys():
                acc_dict[acc] = 1
                hash_start= False
            else:
                acc_dict[acc] = acc_dict[acc] + 1
                hash_start = False
        #Swap these hash starts to use all of the records in the file
        #elif hash_start == True and line[0] == "#":
        #    hash_start = False
        #    continue

    acc_dict_sorted = sorted(
        acc_dict.items(), key=lambda x: x[1], reverse=True)
    return acc_dict_sorted


def initiate_database():
    column_names = ['Accession', 'TaxID', 'Scientific']
    database = pd.DataFrame(columns=column_names)
    return database


def query_ncbi(accesion_list, current_database):
    sci_dict = {}
    count = 0

    unzipped = [[i for i, j in accesion_list],
                [j for i, j in accesion_list]]

    pd.DataFrame(unzipped[0]).to_csv("accession.csv",index=False)
    pd.DataFrame(unzipped[1]).to_csv("coint.csv",index=False)
    #print(unzipped[0][0] +' ,'+ ' ,'.join(unzipped[0][1:]))

    tax_handle = Entrez.efetch(db="Nucleotide", id=unzipped[0][0] +' ,'+ ' ,'.join(unzipped[0][1:]), retmode="xml")
    record = Entrez.read(tax_handle)

    tax_id_handle = Entrez.esummary(db="Nucleotide", id=unzipped[0][0] +' ,'+ ' ,'.join(unzipped[0][1:]), retmode="xml")
    tax_records = Entrez.read(tax_id_handle)

    for i in range(0,len(record)): 
        tax_name = record[i]['GBSeq_source']
        tax_id = int(tax_records[i]["TaxId"])
        if tax_name in sci_dict:
            sci_dict[tax_name] = sci_dict[tax_name] + unzipped[1][i]
        else:
            sci_dict[tax_name] = unzipped[1][i]
        new_row = [unzipped[0][i],str(tax_id),tax_name]
        current_database.loc[len(current_database)]= new_row

    return(sci_dict, current_database)
    


    # for item in accesion_list:
    #    acc = item[0]
    #    tax_handle = Entrez.esummary(db="Nucleotide", id=acc, retmode="xml")
    #    tax_records = Entrez.parse(tax_handle)
    #    for tax in tax_records:
    #        tax_id = int(tax["TaxId"])
    #        time.sleep(2/100)
    #        name_handle = Entrez.esummary(db="Taxonomy", id=tax_id, retmode="xml")
    #        name_records = Entrez.parse(name_handle)
    #        for name in name_records:
    #            curr_name = name['ScientificName']
    #            if curr_name not in sci_dict:
    #                sci_dict[curr_name] = item[1]
    #            else:
    #                sci_dict[curr_name] = sci_dict[curr_name] + item[1]
    #            time.sleep(2/100)
    #        name_handle.close()
    #    tax_handle.close()
    #   count = count + 1
    # Prints how many records have been searched
    # print(count)

    # TODO rebuild using tax_handle = Entrez.efetch(db="Nucleotide", id=acc_list[0] +' ,'+ ' ,'.join(acc_list[1:]), retmode="xml")
    # record = Entrez.read(tax_handle)
    # record[0]['GBSeq_source']

    


def create_batches(accesion_list, batch_size):
    batch_list = []
    for i in range(0, 200):
        first_index = i*batch_size
        second_index = first_index + batch_size
        batch_list.append(accesion_list[first_index:second_index])
    return batch_list


def calc_batch_size(accesion_list):
    size = round(len(accesion_list)/200)
    return size


def main(arg):
    Entrez.email = os.environ.get("EMAIL_KEY")
    Entrez.api_key = os.environ.get("SECRET_KEY")

    blast_name = arg.blast_output
    blast_dict = collect_first_records(blast_name)

    batch_size = calc_batch_size(blast_dict)

    batch_list = create_batches(blast_dict, batch_size)
    
    current_db = initiate_database()

    result_table = {}
    bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
    i=0
    for item in batch_list:
       
        bar.update(i)
        res_table, current_db = query_ncbi(item, current_db)
        for item in res_table.keys():
            if item in result_table:
                result_table[item] = result_table[item] + res_table.get(item)
            else:
                result_table[item] = res_table.get(item)
        #print(result_table) 
        i = i+1
        
    print(result_table)
    print(current_db)
    print(sum(result_table.values()))

    # print(batch_list)

    #sci_name_dict = query_ncbi(blast_dict, empty_db)


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


import xml.etree.ElementTree as ET
sci_req = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&query_key=&id=P023082,AC279172,KT184333,KT184354,MK810783,CP034425,MN307743&retmode=xml').content
tax_req = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&query_key=&id=P023082,AC279172,KT184333,KT184354,MK810783,CP034425,MN307743&retmode=xml').text
sci_tree = ET.ElementTree(ET.fromstring(sci_req))
tax_tree = ET.ElementTree(ET.fromstring(tax_req))
sci_root = sci_tree.getroot()
tax_root = tax_tree.getroot()
#will get the species names 
for child in sci_root:
    child.findtext('GBSeq_organism')

soup = BeautifulSoup(tax_req, 'lxml')
res = soup.findAll("item")
for i in range(1, len(res)):
    if res[i]['name'] == "TaxId":
       curr_id = res[i].text

