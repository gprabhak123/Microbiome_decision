#!/usr/bin/env python

import csv
import numpy as np
import pandas as pd
import argparse
import pickle as pic
import os
import math
import subprocess
#from Bio.SearchIO import FastaIO
from Bio import SearchIO

master_table_columns = ['Hit_tcid', 'Hit_xid', 'CE', 'hit_tms_no', 'Genome1', 'Genome2']
template = pd.DataFrame(columns=master_table_columns)
sw_output = {}
genomes = []
matches = {}
query_tms = {}
query_to_run = []

tms_info = {}

def get_matches(directory):
    files = os.listdir(directory + '/greens')
    
    for file in files:
        df = pd.read_table(os.path.join(directory + '/greens', file))
        
        genome = file.split('.')[0] + '.' + file.split('.')[1]
        genome = genome.split('_')[0] + '_' + genome.split('_')[1]
        matches[genome] = {}
        
        protein_tcid_dict = {}

        with open(directory + '/' + genome + '/tcdbprots.txt', 'r') as f:
            lines = f.readlines()

        for line in lines:
            tc_acc = line.strip()
            protein_tcid_dict[tc_acc.split('-')[0]] = tc_acc.split('-')[1]

        for index, row  in df.iterrows():
            tc_acc = ''
            if row['Hit_tcid'] not in protein_tcid_dict:
                xid = row['Hit_xid']
                fix_cmd = f'grep {xid} ~/db/blastdb/tcdb.faa'
                correct_tcdb = subprocess.run([fix_cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True).stdout.strip().split('>')[1]
                tc_acc = correct_tcdb
            else:
                tc_acc = row['Hit_tcid'] + '-' + row['Hit_xid']


            if tc_acc in matches[genome]:
                matches[genome][tc_acc].append(row['#Query_id'])
            else:
                matches[genome][tc_acc] = [row['#Query_id']]
            

            #query_tms[row['#Query_id']] = row['Query_n_TMS']
            #print(type(row['Query_n_TMS']))

def is_float(string):
    try:
        float(string)
        return True

    except ValueError:
        return False

'''
def extract_pickle(hmmtop_file, genome):
    query_data = {}
    target_data = {}

    with open(hmmtop_file, "rb") as file:
        # Deserialize the data using pickle.load()
        data = pic.load(file)

    # Extracts all query data and target data from hmmtop file
    query_data = data['queries']

    for key in data['tcdb']:
        key_elems = key.split('|')
        new_key = key_elems[3] + '-' + key_elems[2].split('.')[0]

        target_data[new_key] = data['tcdb'][key]
    
    tms_info[genome] = {'target_data': target_data, 'query_data': query_data}
'''
# PRE-CONDITION: all green files are located within a directory titled greens/ and there is an hmmtop.db file for each genome being worked with
def getSmithWaterman(directory):
    row_to_add = template.copy()
    data_to_add = {}

    files = os.listdir(directory + '/greens')

    for file in files:
        genome = file.split('.')[0] + '.' + file.split('.')[1]
        genome = genome.split('_')[0] + '_' + genome.split('_')[1] 
        genomes.append(genome)
        

        if os.path.exists(directory + '/' + genome) and os.path.isdir(directory + '/' + genome):
            continue
        tcdb_proteins = []
        query_proteins = []
        df = pd.read_table(os.path.join(directory + '/greens', file))
    
        gz_cmd = f'mkdir {directory}/{genome}; gunzip -c /ResearchData/Microbiome/Assemblies/{genome}/*.faa.gz > {directory}/{genome}/{genome}.faa'


        print(gz_cmd)
        os.system(gz_cmd)
        
        new_header = True
        with open(directory + '/' + genome + '/' + genome + '.faa', 'r') as r:
            with open(directory + '/' + genome + '/' + genome + '_edited.faa', 'w') as w:
                for line in r:
                    if '>' in line:
                        new_header = True

                    if '>lcl' in line:
                        new_header = False
                    elif new_header == True:
                        w.write(line)

        data_to_add[file] = {}

        for index, row in df.iterrows():
            tcid = row['Hit_tcid']
            protein = row['Hit_xid']

            substrate = row['Predicted_Substrate']
        

            hit_tms = row['Hit_n_TMS']
            query = row['#Query_id']

            tcdb_proteins.append(tcid + '-' + protein)
            query_proteins.append(query)

        with open(directory + '/' + genome + '/tcdbprots.txt', 'w') as f:
            for protein in tcdb_proteins:
                f.write(protein + '\n')
    
        with open(directory + '/' + genome + '/queryprots.txt', 'w') as f:
            for protein in query_proteins:
                if is_float(protein) == False and protein != 'nan' and protein != 'NA':
                    f.write(protein + '\n')
                else:
                    f.write('none' + '\n')

        cmd1 = f'fasta_grepper.py -f  {directory}/{genome}/queryprots.txt  {directory}/{genome}/{genome}_edited.faa > {directory}/{genome}/query.faa'
        cmd2 = f'getseqs tcdb {directory}/{genome}/tcdbprots.txt > {directory}/{genome}/{genome}_tcdb.faa'

        print(cmd1)
        os.system(cmd1)

        print(cmd2)
    
        result = subprocess.run([cmd2], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
        print(result.stderr)
        if len(result.stderr) > 1:
            proteins_fixed = []
            new_additions = {}
            for part in result.stderr.split(' '):
                if '\n' in part:
                    protein = part.split('\n')[0].split('-')[1]
                    if protein in proteins_fixed:
                        continue
                    fix_cmd = f'grep {protein} ~/db/blastdb/tcdb.faa'
                    print(fix_cmd)
                    correct_tcdb = subprocess.run([fix_cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
                    proteins_fixed.append(protein)
                    new_additions[protein] = correct_tcdb.stdout.strip().split('>')[1] 

            print('adding new proteins')
            print(new_additions)
            with open(directory + '/' + genome + '/tcdbprots.txt', 'w') as f:
                for protein in tcdb_proteins:
                    if protein.split('-')[1] in new_additions:
                        f.write(new_additions[protein.split('-')[1]] + '\n')
                    else:
                        f.write(protein + '\n')
        
            print(cmd2)
            os.system(cmd2)


        # smith-waterman call
        sw_cmd = f'ssearch36 -z 11 -k 1000 -W 0 -m 10 {directory}/{genome}/{genome}_tcdb.faa {directory}/{genome}/query.faa > {directory}/{genome}/ssearch.out'
        print(sw_cmd)
        os.system(sw_cmd)
        print(genome + ' done')


# PRE-CONDITION: all green files are labeled with the genome and all genome files have ssearch.out
def parse_sw(directory):
    get_matches(directory)
    dirs = os.listdir(directory)
    matches_recorded = {}
    for d in dirs:
        if os.path.exists(directory + '/' + d) and os.path.isdir(directory + '/' + d):
            if 'ssearch.out' in os.listdir(directory + '/' + d):
                print(d)
                file_path = os.path.join(directory + '/' + d, 'ssearch.out')
                search_results = SearchIO.parse(file_path, "fasta-m10")
                sw_output[d] = {}
                been_matched = []

                matches_recorded[d] = {}
                for search_result in search_results:

                    for query in search_result:
                        #print(query)
                        #exit()
                        #tcacc_len = search_result.seq_len
                        #temp['query_len'] = tcacc_len
                        #temp['query_id'] = search_result.id
                        # QueryResult
                       

                        if search_result.id in sw_output[d] and query.id not in matches[d][search_result.id]:
                            continue
                            
                        comps_found = []
                        for hit in query:
                            # HSP
                            if search_result.id == '3.A.1.134.11-Q9RL74' and d == 'GCF_000013425.1':
                                print('here')

                            temp = {}
                            temp['query_len'] = search_result.seq_len
                            temp['query_id'] = search_result.id

                            temp['hit_len'] = query.seq_len
                            
                            temp['hit_id'] = hit.hit.id
                            temp['eval'] = hit.evalue
                            temp['pident'] = int(hit.ident_pct)
                            temp['query_seq'] = str(hit.query.seq)
                            temp['query_start'] = hit.query_start + 1
                            temp['query_end'] = hit.query_end
                            temp['hit_seq'] = str(hit.hit.seq)
                            temp['hit_start'] = hit.hit_start + 1
                            temp['hit_end'] = hit.hit_end
                            


                            if search_result.id in sw_output[d] and hit.hit.id in matches[d][search_result.id] and hit.hit.id not in comps_found:
                                for match in matches[d][search_result.id]:
                                    if match == hit.hit.id and hit.hit.id not in matches_recorded[d][search_result.id] and temp['hit_id'] not in been_matched and float(temp['eval']) < 1e-3:
                                        sw_output[d][search_result.id].append(temp)
                                        matches_recorded[d][search_result.id].append(temp['hit_id'])
                                        been_matched.append(temp['hit_id'])
                                        comps_found.append(temp['hit_id'])

                            if search_result.id not in sw_output[d] and float(temp['eval']) < 1e-3 and hit.hit.id in matches[d][search_result.id]:
                                sw_output[d][search_result.id] = [temp]
                                matches_recorded[d][search_result.id] = [temp['hit_id']]
                                been_matched.append(temp['hit_id'])
                                comps_found.append(temp['hit_id'])
                    

'''
# TODO: ask R2 if we can still use the info from hmmtop.db as that was from blast results
def getTMOverlap():
    for genome in genomes:
        sw_out = sw_output[genome]
        
        mmseqs = {}
        target_data = tms_info[genome]
        query_data = tms_info[genome]
        hmmtop_dict = {}
        for query in sw_out:
            if query not in mmseqs:
                mmseqs[query] = {'qaln': sw_out[query]['query_seq'], 'taln': sw_out[query]['hit_seq'], 'target': sw_out[query]['hit_id'], 'qstart': sw_out[query]['query_start'], 'qend': sw_out[query]['query_end'], 'tstart': sw_out[query]['hit_start'], 'tend': sw_out[query]['hit_end']}

                qtms = {}
                if query in query_data:
                    qtms['tms'] = list(query_data[qid].values())


                    tms_dict[query] = len(qtms['tms'])
                    if query not in hmmtop_dict:
                        hmmtop_dict[query] = qtms
'''

# Given a query it returns its data in df format
def getInfoAsRow(genome, query):
    #TODO: ask Rif if he wants data as a dict or df
    
    if genome not in sw_output:
        print('Genome not been processed: ' + genome)
        return

    genome_data = sw_output[genome]
    
    

    if query not in genome_data:
        raise Exception('Query was not processed: ' + query)
        return 
    
    
    query_data = genome_data[query]

    df = pd.DataFrame(query_data)

    return df



#getSmithWaterman('/Users/gautham/microbiome_project/test_comparisons')
#parse_sw('/Users/gautham/microbiome_project/test_comparisons')
#print(getInfoAsRow('GCF_000013425.1', '3.A.1.134.11-Q9RL74'))


