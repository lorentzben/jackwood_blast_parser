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
import requests
import lxml
from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET
import gc
import pickle
import csv


def collect_first_records(blast_frame):
    blast_file = open(blast_frame, 'r')
    Lines = blast_file.readlines()

    hash_start = False
    acc_dict = {}
    per_q_dict = {}
    per_i_dict = {}
    count=0
    for line in Lines:
        if line[0] == "#":
            hash_start = True
        if hash_start == True and line[0] != "#":
            acc = str.split(line, ",")[1]
            perc_query = str.split(line,",")[5]
            perc_ident = str.split(line,",")[6]
            if acc not in acc_dict.keys():
                acc_dict[acc] = 1
                hash_start= False
                count= count+1
            if acc in acc_dict.keys():
                acc_dict[acc] = acc_dict[acc] + 1
                hash_start = False
                count = count+1
            if acc not in per_q_dict.keys():
                per_q_dict[acc] = list([float(perc_query)])
                hash_start = False
            if acc in per_q_dict.keys():
                per_q_dict[acc].append(float(perc_query))
                hash_start = False
            if acc not in per_i_dict.keys():
                per_i_dict[acc] = list([float(perc_ident)])
                hash_start = False
            if acc in per_i_dict.keys():
                 per_i_dict[acc].append(float(perc_query))
                 hash_start = False
        
        
                 
        #Swap these hash starts to use all of the records in the file
        #elif hash_start == True and line[0] == "#":
        #    hash_start = False
        #    continue

    gc.collect()
    print("record count: "+str(count))
    acc_dict_sorted = sorted(
        acc_dict.items(), key=lambda x: x[1], reverse=True)
    return [acc_dict_sorted, per_q_dict, per_i_dict]


def initiate_database():
    column_names = ['Accession', 'TaxID', 'Scientific']
    database = pd.DataFrame(columns=column_names)
    return database



def query_ncbi_xml(accesion_list, current_database, api_key):
    db = current_database
    sci_dict = {}
    to_be_removed = []
    acc_dict = dict(accesion_list)
    for item in accesion_list:
        curr_acc = item[0]
        curr_value = item[1]
        
        #print(curr_acc in current_database["Accession"].values)
        if curr_acc in current_database["Accession"].values:
            #print(curr_acc +" found in db")
            curr_name = str(current_database[current_database["Accession"].values == curr_acc]["Scientific"].values)
            to_be_removed.append(curr_acc)
            if curr_name not in sci_dict:
                sci_dict[curr_name] = curr_value
            else:
                sci_dict[curr_name] = sci_dict[curr_name] + curr_value
    
    for item in to_be_removed:
        del acc_dict[item] 

    accesion_list = list(acc_dict.items())
    #print(accesion_list)
    unzipped = [[i for i, j in accesion_list],
                [j for i, j in accesion_list]]
    #print(unzipped)
    if unzipped != [[],[]]:
        
        try:
            id_list = unzipped[0][0] +','+ ','.join(unzipped[0][1:])
        except IndexError:
            try:
                id_list = unzipped[0][0]
            except IndexError:
                pass

        #print(id_list)
        gc.collect()
        
        try:
            sci_req = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&query_key=&id='+id_list+'&retmode=xml').content
        except ChunkEncodingError as ex:
            print(f"Invalid chunk encoding {str(ex)}")
            exit(1)
        try:
            tax_req = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&query_key=&id='+id_list+'&retmode=xml').content
        except ChunkEncodingError as ex2:
            print(f"Invalid chunk encoding {str(ex)}")
            exit(1)
        
    

        count = 0
        sci_names = []
        soup_sci = BeautifulSoup(sci_req,'lxml')
        sci_res = soup_sci.findAll("gbseq_organism")
        for item in sci_res:
            curr_name = item.text
            sci_names.append(curr_name)
            if curr_name not in sci_dict:
                sci_dict[curr_name] = unzipped[1][count]
            else:
                sci_dict[curr_name] = sci_dict[curr_name] + unzipped[1][count]
            count = count + 1
        gc.collect()
            
        soup_tax = BeautifulSoup(tax_req, 'lxml')
        tax_res = soup_tax.findAll("item")
        i = 0
        j = 0
        for item in tax_res:
            if tax_res[i]['name'] == "TaxId":
                #print("i, j: " +str(i)+" "+str(j))
                curr_id = tax_res[i].text
                #print("current_id: "+curr_id)
                new_row = pd.DataFrame({'Accession':[unzipped[0][j]], 
                'TaxID':[str(curr_id)], 
                'Scientific':[sci_names[j]]})
                #print(new_row)
                db = db.append(new_row, ignore_index=True)
                j = j+1
            i = i+1
        gc.collect()
    return(sci_dict, db)
    

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
    

def create_batches(accesion_list, batch_size, divider):
    batch_list = []
    for i in range(0, divider):
        first_index = i*batch_size
        second_index = first_index + batch_size
        batch_list.append(accesion_list[first_index:second_index])
    gc.collect()
    return batch_list


def calc_batch_size(accesion_list):
    print("accession_count: "+str(len(accesion_list)))
    divider = 200
    size = round(len(accesion_list)/divider)
    return size, divider

def write_db_to_file(current_database,db_name):
    current_database.to_csv(db_name,index=False)

def write_res_tab_to_file(res_table,res_name):
    names = list(res_table.keys())
    values = list(res_table.values())
    name_list = []
    for item in names:
        fix_name = item.strip('][').split(', ')[0].replace("'",'')
        name_list.append(fix_name)
    
    res_df = pd.DataFrame([name_list, list(res_table.values())]).T
    res_df.to_csv(res_name,index=False)
    return res_df

def g_mean(x):
    # from https://www.statology.org/geometric-mean-python/
    #input: list of int != 0
    #output: geometric mean of list
    x =[i for i in x if i !=0]
    a = np.log(x)
    return round(np.exp(a.mean()),2)


def convert_acc_query_coverage_dict_to_sci_qc_dict(database, query_acc_dict):
    database = database.sort_values(by='Scientific')
    sci_qc_dict = {}
    missing_acc = []

    #qafile = open('query_acc_dict.bjl', 'ab')
    #pickle.dump(query_acc_dict,qafile)
    #qafile.close()

    for item in database.iterrows():
        sci_name = item[1][2]
        acc = item[1][0]
        if sci_name not in sci_qc_dict:
            if query_acc_dict.get(acc) != None:
                sci_qc_dict[sci_name] = query_acc_dict.get(acc)
        elif sci_name in sci_qc_dict:
            try:
                val_list = query_acc_dict.get(acc)
                #print("valist: " + str(val_list))
                for value in val_list:
                    #print("value: " + str(value))
                    list(sci_qc_dict[sci_name]).append(value)
            except TypeError as ex3:
                missing_acc.append(acc)

    #mafile = open('missing_acc.bjl', 'ab')
    #pickle.dump(missing_acc, mafile)
    #mafile.close()

    #dbfile = open("database.bjl", 'ab')
    #pickle.dump(database,dbfile )
    #dbfile.close()
    return sci_qc_dict, missing_acc

def convert_acc_identify_dict_to_sci_ident_dict(database, ident_acc_dict):
    database = database.sort_values(by='Scientific')
    sci_id_dict = {}
    missing_acc = []

    #idfile = open('id_acc_dict.bjl', 'ab')
    #pickle.dump(query_acc_dict,idfile)
    #idfile.close()

    for item in database.iterrows():
        sci_name = item[1][2]
        acc = item[1][0]
        if sci_name not in sci_id_dict:
            if ident_acc_dict.get(acc) != None:
                sci_id_dict[sci_name] = ident_acc_dict.get(acc)
        elif sci_name in sci_id_dict:
            try:
                val_list = ident_acc_dict.get(acc)
                #print("valist: " + str(val_list))
                for value in val_list:
                    #print("value: " + str(value))
                    list(sci_id_dict[sci_name]).append(value)
            except TypeError as ex3:
                missing_acc.append(acc)
    #mafile = open('missing_acc.bjl', 'ab')
    #pickle.dump(missing_acc, mafile)
    #mafile.close()

    #dbfile = open("database.bjl", 'ab')
    #pickle.dump(database,dbfile )
    #dbfile.close()
    return sci_id_dict, missing_acc


def calc_average_query_coverage_per_species (query_sci_dict):
    average_q_dict = {}
    for item in query_sci_dict.keys():
        query_list = query_sci_dict[item]
        average_q_dict[item] = g_mean(query_list)
    return(average_q_dict)
        

def calc_average_identity_coverage_per_species (id_sci_dict):
    average_id_dict = {}
    for item in id_sci_dict.keys():
        id_list = id_sci_dict[item]
        average_id_dict[item] = g_mean(id_list)
    return(average_id_dict)

def main(arg):
    Entrez.email = os.environ.get("EMAIL_KEY")
    Entrez.api_key = os.environ.get("SECRET_KEY")

    blast_name = arg.blast_output
    blast_dict, query_coverage_dict, identity_dict = collect_first_records(blast_name)

    #print(query_coverage_dict)
    #print(identity_dict)

    batch_size , divider = calc_batch_size(blast_dict)

    batch_list = create_batches(blast_dict, batch_size, divider)
    
    if arg.database:
        current_db = pd.read_csv(arg.database)
    else:
        current_db = initiate_database()
    
    if arg.database_file:
        db_name = arg.database_file
    else:
        db_name = "database_1.csv"

    if arg.result_file:
        res_name = arg.result_file
    else:
        res_name = "result_1.csv"
        

    result_table = {}
    bar = progressbar.ProgressBar(max_value=divider)
    i=0
    for item in batch_list:
        bar.update(i)
        res_table, current_db = query_ncbi_xml(item, current_db, Entrez.api_key)
        for item in res_table.keys():
            if item in result_table:
                result_table[item] = result_table[item] + res_table.get(item)
            else:
                result_table[item] = res_table.get(item)
        #print(result_table) 
        #print(current_db)
        gc.collect()
        i = i+1
        
    #print(result_table)
    #print(current_db)
    #print(sum(result_table.values()))

    write_db_to_file(current_db,db_name)
    res_df = write_res_tab_to_file(result_table,res_name)

    sci_query_coverage_dict, missing_qc_acc = convert_acc_query_coverage_dict_to_sci_qc_dict(current_db, query_coverage_dict)
    
    sci_id_coverage_dict, missing_id_acc = convert_acc_identify_dict_to_sci_ident_dict(current_db, identity_dict)

    
    #with open('test.csv', 'w') as f:
    #    for key in sci_query_coverage_dict.keys():
    #        f.write("%s, %s\n" % (key, sci_query_coverage_dict[key]))

    query_coverage_dict = calc_average_query_coverage_per_species(sci_query_coverage_dict)
    identity_dict = calc_average_identity_coverage_per_species(sci_id_coverage_dict)
    
    
    qc_id_table = pd.concat([pd.DataFrame(query_coverage_dict, index=[0]).T,pd.DataFrame(identity_dict, index=[0]).T], axis=1)
    qc_id_table.columns = ['ave_query_coverage','ave_identitiy']

    res_df.columns = ['taxa','count']
    res_df.set_index('taxa', inplace=True)
    qc_id_table.index.name = 'taxa'

    joined_table = res_df.join(qc_id_table, how='outer')

    joined_table.to_csv(res_name[:-4]+"_more_info.csv",index=True)

   
    # print(batch_list)

    


if __name__ == "__main__":
    # Build Argument Parser in order to facilitate ease of use for user
    parser = argparse.ArgumentParser(
        description="Digests blast format 6 outputs and reports counts of each taxa identified")
    parser.add_argument('-b', '--blast_output', action='store', required=True,
                        help="name of the blast output file", dest='blast_output')
    parser.add_argument('-d', '--db', action='store', required=False,
                        help="name of the database from previous runs", dest='database')
    parser.add_argument('-o','--dbo', action='store',required=False, help="name of the database to save to disk", dest='database_file')
    parser.add_argument('-n','--res_o', action='store',required=False, help="name of the result to save to disk", dest='result_file')

    args = parser.parse_args()
    main(args)


