#!/usr/bin/env python3

import Bio
from lxml import objectify as xml_objectify
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
from collections import defaultdict
import warnings

warnings.filterwarnings(action='ignore', category=FutureWarning)



def xml_to_dict(xml_str):
    """ Convert xml to dict, using lxml v3.4.2 xml processing library from user radtek on Stack Overflow """
    def xml_to_dict_recursion(xml_object):
        dict_object = xml_object.__dict__
        if not dict_object:
            return xml_object
        for key, value in dict_object.items():
            dict_object[key] = xml_to_dict_recursion(value)
        return dict_object
    return xml_to_dict_recursion(xml_objectify.fromstring(xml_str))


def retrieve_strain_from_list_of_dict(xml_dict_list):
    # Takes in a list of dicts from gbqualifier
    strain_list = []
    for dict in xml_dict_list:
        try:
            curr_key = str(list(dict.values())[0])
            curr_val = str(list(dict.values())[1])
        except IndexError:
            curr_key = ""
            curr_val = ""
            pass
        if (curr_key in ["strain","isolate"]) :
            strain_list.append(curr_val)
    return strain_list

def update_strain_tuples(strain, count, strain_count_list):
    updated_list_of_tuples = []

    unzipped = [[i for i, j in strain_count_list],
                [j for i, j in strain_count_list]]

    if strain in unzipped[0]:
        loc = unzipped[0].index(strain)
        unzipped[1][loc] = unzipped[1][loc] + count
    else:
        unzipped[0].append(strain)
        unzipped[1].append(count)
    
    updated_list = list(zip(unzipped[0],unzipped[1]))
    return updated_list


def collect_first_records(blast_frame):
    # Loads the blast result file to memory
    blast_file = open(blast_frame, 'r')
    line = blast_file.readlines()

    acc_dict = {}
    per_q_dict = {}
    per_i_dict = {}
    count=0

    record_list_of_list = []
    chunk = []
    last_hash = True
    for record in line:
        
        #print("last_hash: ", last_hash)
        if not record.startswith('#') and last_hash :
            chunk.append(record)
            last_hash = False
        if not record.startswith('#') and not last_hash:
            last_hash = False
        if record.startswith('#') and last_hash:
            last_hash = True
        if record.startswith('#') and not last_hash:
            record_list_of_list.append(chunk)
            chunk = []
            last_hash = True

    '''
    # This method code is used if you want all records to be examined, it will make you results inflated approx 10x based on blast hits.
    last_hash = True
    chunk = []
    for record in line:
        
        #print("last_hash: ", last_hash)
        if not record.startswith('#') and last_hash :
            last_hash = False
        if not record.startswith('#') and not last_hash:
            chunk.append(record)
            last_hash = False
        if record.startswith('#') and last_hash:
            last_hash = True
        if record.startswith('#') and not last_hash:
            record_list_of_list.append(chunk)
            chunk = []
            last_hash = True
    '''

    # instantiates dicts for accessions, and coverage stats
    acc_dict = {}
    per_q_dict = {}
    per_i_dict = {}
    count=0
    # iterates over each line and uses the hash from comments to determine which record is first
    for record in record_list_of_list:
        #print(record)
        for line in record:
            # since this is a csv like document we can select indicies as values
            acc = str.split(line, ",")[1]
            perc_query = str.split(line,",")[5]
            perc_ident = str.split(line,",")[6]
            #print("acc: ", acc," perc_query: ", perc_query, " perc_ident: ", perc_ident)
            # checks if value exists in the accession dict and adds it if not 
            if acc not in acc_dict.keys():
                acc_dict[acc] = 1
                count= count + 1
                #print('added accession: ', count)
            # checks if value exists in the accession dict and increases count 
            elif acc in acc_dict.keys():
                acc_dict[acc] = acc_dict[acc] + 1
                count = count + 1
                #print('updated accession: ', count)
            # performs same check as accessions but creates a list of coverages
            if acc not in per_q_dict.keys():
                per_q_dict[acc] = list([float(perc_query)])
            elif acc in per_q_dict.keys():
                per_q_dict[acc].append(float(perc_query))
            if acc not in per_i_dict.keys():
                per_i_dict[acc] = list([float(perc_ident)])
            elif acc in per_i_dict.keys():
                per_i_dict[acc].append(float(perc_query))
       
    
    gc.collect()
    # sorts accession dict so that the most popular accessions are at top to speed up quering.
    # if the record count does not match # BLAST processed XXXX queries it is because you had queries with 0 hits, do not worry
    print("record count: ", count)
    acc_dict_sorted = sorted(
        acc_dict.items(), key=lambda x: x[1], reverse=True)
    #print(acc_dict_sorted)
    return [acc_dict_sorted, per_q_dict, per_i_dict]

def initiate_database():
    # intiates an empty dataframe with headers
    column_names = ['Accession', 'TaxID', 'Scientific','Strain']
    database = pd.DataFrame(columns=column_names)
    return database



def query_ncbi_xml(accesion_list, current_database, api_key,verb):
    db = current_database
    sci_dict = {}
    strain_dict = {}
    to_be_removed = []
    acc_dict = dict(accesion_list)
    # checks if an api key was provided and gives a dummy value if not included
    if (api_key is None) or (api_key == ""):
        api_key = ""
    else:
        api_key = "&api_key="+api_key
    # iterates over each accession and count and separates them out
    for item in accesion_list:
        curr_acc = item[0]
        curr_value = item[1]
        
        #print(curr_acc in current_database["Accession"].values)
        if curr_acc in current_database["Accession"].values:
            
            # selects current scientific name from the database if it exists
            curr_name = " ".join(str(x) for x in list(current_database[current_database["Accession"].values == curr_acc]["Scientific"].values))
            curr_strain = " ".join(str(y) for y in list(current_database[current_database["Accession"].values == curr_acc]["Strain"].values))
            
            # add the accession number to a list to be removed from the batch list
            to_be_removed.append(curr_acc)
            # check if the current scientific name is in the dictionary and add or update value
            if curr_name not in sci_dict:
                
                sci_dict[curr_name] = [(curr_strain,curr_value)]
            else:
                
                updated_res = update_strain_tuples(curr_strain,curr_value, sci_dict.get(curr_name))
                
                sci_dict[curr_name] = update_strain_tuples(curr_strain,curr_value, sci_dict.get(curr_name))
            
            if str(curr_strain) not in ['["Unkown"]',"['Unkown']"]:
                if curr_strain not in strain_dict:
                    strain_dict[curr_strain] = curr_value
                else:
                    strain_dict[curr_strain] = strain_dict[curr_strain] + curr_value  

    # remove the accession from the accession list to be queried if its not in the database. 
    for item in to_be_removed:
        del acc_dict[item] 

    # get a list of the accessions to be queried
    accesion_list = list(acc_dict.items())
    #print(accesion_list)
    # separates list of accessions and list of counts into a nested list 
    unzipped = [[i for i, j in accesion_list],
                [j for i, j in accesion_list]]
    #print(unzipped)
    # as long as there is something in both lists do this 
    if unzipped != [[],[]]:
        
        # this is in place in case there is only one accession ID passed through 
        # first case joins all accessions with commas in between
        # second takes one acession (odd number)
        # final is if there are no more accessions to let the program continue gracefully (even number)
        try:
            id_list = unzipped[0][0] +','+ ','.join(unzipped[0][1:])
        except IndexError:
            try:
                id_list = unzipped[0][0]
            except IndexError:
                pass
        if(verb):
            print(id_list)
        gc.collect()
        # query for scientific names from list of accessions
        if(verb):
            print('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore'+api_key+'&id='+id_list+'&retmode=xml')
        try:
            sci_req = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore'+api_key+'&id='+id_list+'&retmode=xml').content
        except ChunkedEncodingError as ex:
            print(f"Invalid chunk encoding {str(ex)}")
            #TODO add these requests to a list and then retry 
            exit(1)
        # query for taxID from list of accessions 
        try:
            tax_req = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore'+api_key+'&id='+id_list+'&retmode=xml').content
        except ChunkedEncodingError as ex2:
            print(f"Invalid chunk encoding {str(ex2)}")
    
            exit(1)
        
    
        count = 0
        sci_names = []

        # from https://stackoverflow.com/a/14463362

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            soup_sci = BeautifulSoup(sci_req,'lxml')

        # Strain isolate information gathering
        qualifiers = soup_sci.findAll("gbqualifier")

        # turns gbqualifier names and values into a dict
        xml_dict = [xml_to_dict(str(elem)) for elem in qualifiers]

        # flattens a list of dicts where the keys are all the same in both dicts
        strainlist = retrieve_strain_from_list_of_dict(xml_dict)
        

        #strainlist index is in the same order as accessions
                
        # Beautiful Soup is a parser to make sense of the XML 
        # this selects for all of the organism names from NCBI so we can link them to accession IDs
        
        sci_res = soup_sci.findAll("gbseq_organism")
        # connects counts of accession to scientific name and adds it to result list or updates the value
        # TODO nest the strain check here like 164 to add or updated tuple values
        p = 0 
        for item in sci_res:
            #print("sci res: " +str(list(sci_res)))
            #print("strain list: " + str(list(strainlist)))
            curr_name = item.text
            sci_names.append(curr_name)
            try:
                curr_strain = strainlist[p]
            except IndexError:
                curr_strain = "Unknown"

            if curr_name not in sci_dict:
                sci_dict[curr_name] = [(curr_strain, unzipped[1][count])] 
            else:
                sci_dict[curr_name] = update_strain_tuples(curr_strain,unzipped[1][count], sci_dict.get(curr_name)) 
            
            if curr_strain not in strain_dict:
                strain_dict[curr_strain] = unzipped[1][count]
            else:
                strain_dict[curr_strain] = strain_dict[curr_strain] + unzipped[1][count] 
            count = count + 1
            p = p+1
        gc.collect()

        #from https://stackoverflow.com/a/14463362
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")    
            soup_tax = BeautifulSoup(tax_req, 'lxml')
        tax_res = soup_tax.findAll("item")
        i = 0
        j = 0
        
        # links tax id to accession number 
        #tax_res[i]
        for item in tax_res:
            #print("item:" , item, " i: ",i)
            if item['name'] == "TaxId":
                #print("i, j: " +str(i)+" "+str(j))
                curr_id = tax_res[i].text
                try:
                    curr_strain = strainlist[j]
                except IndexError:
                    curr_strain = "Unkown"
                #print("current_id: "+curr_id)
                #print("length scinames", len(sci_names))
                #print("length strain", len(strainlist))

                new_row = pd.DataFrame({'Accession':[unzipped[0][j]], 
                'TaxID':[str(curr_id)], 
                'Scientific':[sci_names[j]],
                'Strain':[curr_strain]})
                db = db.append(new_row, ignore_index=True)
                j = j+1
                #print(new_row)
                
            i = i+1
        gc.collect()
        #print(sci_dict)
        #print(db)
    
        #print(strain_dict)
    return(sci_dict, db, strain_dict)
    

def query_ncbi(accesion_list, current_database):
    # TODO Not actually implemented, was a intermediary development method
    sci_dict = {}
    count = 0

    # splits the accession number and count
    unzipped = [[i for i, j in accesion_list],
                [j for i, j in accesion_list]]

    # saves accession and count to disk (can be turned off)
    pd.DataFrame(unzipped[0]).to_csv("accession.csv",index=False)
    pd.DataFrame(unzipped[1]).to_csv("count.csv",index=False)
    #print(unzipped[0][0] +' ,'+ ' ,'.join(unzipped[0][1:]))

    # this reads the query into memory similar to how you open a file in python for scientific names
    tax_handle = Entrez.efetch(db="Nucleotide", id=unzipped[0][0] +' ,'+ ' ,'.join(unzipped[0][1:]), retmode="xml")
    record = Entrez.read(tax_handle)

    # this reads the query into memory similar to how you open a file in python for tax ids
    tax_id_handle = Entrez.esummary(db="Nucleotide", id=unzipped[0][0] +' ,'+ ' ,'.join(unzipped[0][1:]), retmode="xml")
    tax_records = Entrez.read(tax_id_handle)

    # goes through accession resutlts and constructs a result tabl 
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

    '''
    if batch_size*divider < len(accesion_list): 
        divider = len(accesion_list)
        size = 1
    '''

    batch_list = []
    # subsets long list by the provided divider and returns a list
    for i in range(0, divider):
        first_index = i*batch_size
        second_index = first_index + batch_size
        #print("first index: ", first_index, " second index ", second_index)
        batch_list.append(accesion_list[first_index:second_index])
    #print("leftover indicies",batch_size*divider, ' : ',len(accesion_list))
    batch_list.append(accesion_list[batch_size*divider:len(accesion_list)])
    #print("accession list: ", accesion_list)
    #print("batch list: ", batch_list)
    gc.collect()
    return batch_list


def calc_batch_size(accesion_list, acc_len):
    # determines how many batches can be constructed based on divider
    print("accession_count: "+str(len(accesion_list)))
    
    divider = 200
    if len(accesion_list) > divider:
        size = round(len(accesion_list)/divider)
    else:
        size = len(accesion_list)
        

    #print("size: ", size, " divider: ", divider)
    return size, divider

def write_db_to_file(current_database,db_name):
    # wrapper function to save database to csv file
    current_database.to_csv(db_name,index=False)

def write_res_tab_to_file(res_table,res_name):
    # turns a dict into a dataframe to save as a csv file
    names = list(res_table.keys())
    values = list(res_table.values())
    

    column_names = ['Scientific','Strain', 'Count']
    res_df = pd.DataFrame(columns=column_names)
   
    new_row = []
    for i in range(0,len(names)):
        curr_name = names[i]
        for j in range(0,len(values[i])):
            curr_strain = values[i][j][0]
            curr_count = values[i][j][1]
            #print(curr_name," ", curr_strain, " ",curr_count)
            new_row = pd.DataFrame({'Scientific': [curr_name],
                'Strain': [curr_strain], 
                'Count': [curr_count]})
            res_df = res_df.append(new_row, ignore_index=True)
           
            
    res_df.to_csv(res_name,index=False)
    return res_df

def dict_list_tuples_to_dataframe(dict_to_mod, third_col):
    # turns a dict into a dataframe to save as a csv file
    names = list(dict_to_mod.keys())
    values = list(dict_to_mod.values())

    column_names = ['taxa','strain', third_col]
    res_df = pd.DataFrame(columns=column_names)

    new_row = []
    for i in range(0,len(names)):
        curr_name = names[i]
        for j in range(0,len(values[i])):
            curr_strain = values[i][j][0]
            curr_count = values[i][j][1]
            #print(curr_name," ", curr_strain, " ",curr_count)
            new_row = pd.DataFrame({'taxa': [curr_name],
                'strain': [curr_strain], 
                third_col:  [curr_count]})
            res_df = res_df.append(new_row, ignore_index=True)
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

    # iterates through the database and checks for accessions, either adds or updates query coverage values

    for item in database.iterrows():
        
        sci_name = item[1][2]
        strain = item[1][3]
        acc = item[1][0]
        #print("acc " , acc, " acc dict get ", query_acc_dict.get(acc))
        if sci_name not in sci_qc_dict:
            if strain != "Unknown":
                if query_acc_dict.get(acc) != None:
                    sci_qc_dict[sci_name] = [(strain, query_acc_dict.get(acc))]
            else:
                if query_acc_dict.get(acc) != None:
                    sci_qc_dict[sci_name] = [(strain, query_acc_dict.get(acc))]
        elif sci_name in sci_qc_dict:
            try:
                if strain != "Unknown":
                    sci_qc_dict[sci_name] = update_strain_tuples(strain, query_acc_dict.get(acc), sci_qc_dict.get(sci_name))
                else:
                    sci_qc_dict[sci_name] = [("Unknown", query_acc_dict.get(acc))]
                
            except TypeError as ex3:
                missing_acc.append(acc)

    #mafile = open('missing_acc.bjl', 'ab')
    #pickle.dump(missing_acc, mafile)
    #mafile.close()

    #dbfile = open("database.bjl", 'ab')
    #pickle.dump(database,dbfile )
    #dbfile.close()
    #print(sci_qc_dict)
    return sci_qc_dict, missing_acc

def convert_acc_identify_dict_to_sci_ident_dict(database, ident_acc_dict):
    database = database.sort_values(by='Scientific')
    sci_id_dict = {}
    missing_acc = []

    #idfile = open('id_acc_dict.bjl', 'ab')
    #pickle.dump(query_acc_dict,idfile)
    #idfile.close()

    # iterates through the database and checks for accessions, either adds or updates identity coverage values

    for item in database.iterrows():
        sci_name = item[1][2]
        strain = item[1][3]
        acc = item[1][0]
        #print("acc ", acc, " ident acc get ",ident_acc_dict.get(acc))
        if sci_name not in sci_id_dict:
            if strain != "Unknown":
                if ident_acc_dict.get(acc) != None:
                    sci_id_dict[sci_name] = [(strain, ident_acc_dict.get(acc))]
            else:
                if ident_acc_dict.get(acc) != None:
                    sci_id_dict[sci_name] = [(strain, ident_acc_dict.get(acc))]
                #sci_id_dict[sci_name] = ident_acc_dict.get(acc)
        elif sci_name in sci_id_dict:
            try:
                if strain != "Unknown":
                    sci_id_dict[sci_name] = update_strain_tuples(strain, ident_acc_dict.get(acc), sci_id_dict.get(sci_name))
                else:
                    sci_id_dict[sci_name] = [("Unknown", ident_acc_dict.get(acc))]
            except TypeError as ex3:
                missing_acc.append(acc)
    #mafile = open('missing_acc.bjl', 'ab')
    #pickle.dump(missing_acc, mafile)
    #mafile.close()

    #dbfile = open("database.bjl", 'ab')
    #pickle.dump(database,dbfile )
    #dbfile.close()
    #print(sci_id_dict)
    return sci_id_dict, missing_acc


def calc_average_query_coverage_per_species (query_sci_dict):
    # iterates over each specie name and calculates the geometric mean for query coverage
    average_q_dict = {}
    for item in query_sci_dict.keys():
        #print(item)
        query_list = query_sci_dict[item]
        for strain in query_list:
            #print("strain: ", strain[0], " value: ", strain[1])
            if item not in average_q_dict:
                #print("adding ", strain[0], " value ", g_mean(strain[1]))
                average_q_dict[item] = [(strain[0], g_mean(strain[1]))]
            elif item in average_q_dict:
                #print("updating ", strain[0], " value ", g_mean(strain[1]))
                average_q_dict[item].append((strain[0], g_mean(strain[1])))
    
    return(average_q_dict)
        

def calc_average_identity_coverage_per_species (id_sci_dict):
    # iterates over each specie name and calculates the geometric mean for identity coverage
    average_id_dict = {}
    for item in id_sci_dict.keys():
        id_list = id_sci_dict[item]
        for strain in id_list:
            if item not in average_id_dict:
                average_id_dict[item] = [(strain[0], g_mean(strain[1]))]
            elif item in average_id_dict:
                average_id_dict[item].append((strain[0], g_mean(strain[1])))

    return(average_id_dict)

def main(arg):
    #If user has an NCBI api key, place here for less restricted queurying
    os.environ['SECRET_KEY'] = ""
    os.environ['EMAIL_KEY'] = ""

    Entrez.email = os.environ['EMAIL_KEY']
    Entrez.api_key = os.environ['SECRET_KEY']
   
    
    # input file from params
    blast_name = arg.blast_output
    # Selects first record from each blast result as well as percent query and identity 
    blast_dict, query_coverage_dict, identity_dict = collect_first_records(blast_name)
    #print(blast_dict)
    #print(len(blast_dict))
    # diagnostic print statments
    #print(query_coverage_dict)
    #print(identity_dict)

    # determines batch size so as not to overload querying 
    batch_size , divider = calc_batch_size(blast_dict, len(blast_dict))

    # divides query list into batches
    batch_list = create_batches(blast_dict, batch_size, divider)
    #print(batch_list)
    
    # checks for existing database if none present creates new empty
    if arg.database:
        current_db = pd.read_csv(arg.database)
    else:
        current_db = initiate_database()
    
    # if no database name provided, database_1.csv used 
    if arg.database_file:
        db_name = arg.database_file
    else:
        db_name = "database_1.csv"

    # if not result name file provided result_1.csv used
    if arg.result_file:
        res_name = arg.result_file
    else:
        res_name = "result_1.csv"
    

    # prints info to terminal if selected  
    if args.verbose:
        verbose = True
    else:
        verbose = False
    
    # iterates through batch list and queries ncbi for records, checks if record currently exists and adds or updates 
    result_table = {}
    strain_table = {}
    bar = progressbar.ProgressBar(max_value=len(batch_list))
    i=0
    for item in batch_list:
        bar.update(i)
        res_table, current_db, strain_dict = query_ncbi_xml(item, current_db, Entrez.api_key, verbose)
        #TODO see if this makes sense after getting to the lower logic
        for item in list(res_table.keys()):
            #print("res_table_key :",item)
            #print("res_table_value: ",list(res_table.values()) )
            if item in result_table:
                result_table[item] = result_table[item] + res_table.get(item)
            else:
                result_table[item] = res_table.get(item)

        for item in strain_dict.keys():
            if item in strain_table:
                strain_table[item] = strain_table[item] + strain_table.get(item)
            else:
                strain_table[item] = strain_dict.get(item)
        #print(result_table) 
        #print(current_db)
        # frees some memory 
        gc.collect()
        i = i+1
        
    #print(result_table)
    #print(current_db)
    #print(sum(result_table.values()))

    # save database to disk
    write_db_to_file(current_db,db_name)
    # save query results to disk
    
    res_df = write_res_tab_to_file(result_table,res_name)

    #TODO update with more params
    #strain_df = write_res_tab_to_file(strain_table,"strain_table.csv")

    # calculates geometric mean of average query coverage for each record
    sci_query_coverage_dict, missing_qc_acc = convert_acc_query_coverage_dict_to_sci_qc_dict(current_db, query_coverage_dict)

    # calculates geometric mean of average identity coverage for each record
    sci_id_coverage_dict, missing_id_acc = convert_acc_identify_dict_to_sci_ident_dict(current_db, identity_dict)

    
    #with open('test.csv', 'w') as f:
    #    for key in sci_query_coverage_dict.keys():
    #        f.write("%s, %s\n" % (key, sci_query_coverage_dict[key]))

    # calculates geometric mean of average query coverage for each species
    query_coverage_dict = calc_average_query_coverage_per_species(sci_query_coverage_dict)

    # calculates geometric mean of average identity coverage for each species
    identity_dict = calc_average_identity_coverage_per_species(sci_id_coverage_dict)
    
    # TODO need to use something like write_res_tab_to_file to make two DFs that can be merged on the species then strain cols
    # The result can be merged to the original data, and then written to file

    qc_df = dict_list_tuples_to_dataframe(query_coverage_dict, 'ave_query_coverage')
    id_df = dict_list_tuples_to_dataframe(identity_dict, 'ave_identitiy')

    ave_df = pd.merge(qc_df, id_df,  how='left', left_on=['taxa','strain'], right_on = ['taxa','strain'])
    
    

    # renames columns and sets index to taxa as opposed to numeric
    res_df.columns = ['taxa', 'strain','count']
    res_df.set_index('taxa', inplace=True)

    more_info_df = ave_df = pd.merge(res_df, ave_df,  how='left', left_on=['taxa','strain'], right_on = ['taxa','strain'])
    
    # saves joined table to disk
    more_info_df.to_csv(res_name[:-4]+"_more_info.csv",index=False)

   
    # print(batch_list)

    
# parse parameters from user submission
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
    parser.add_argument('-v', '--verb', action='store_true', required=False, help="verbose output to the terminal", dest='verbose')

    args = parser.parse_args()
    main(args)


