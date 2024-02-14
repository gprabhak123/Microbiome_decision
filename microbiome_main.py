#!/usr/bin/env python

import argparse
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import csv
import numpy as np
import os
import re
import pickle as pic
import subprocess
import subprocess
# from fusion_distribution import isFusion,geneFusions
# from fusion_dist_dir import isFusion,genDict
from fusion import isFusion, geneFusions, genDict, setGenome, check_overlap, fusion_TMS_count
from parseXML import parse
from pprint import pprint
from find_protein import find, extract_all, find_mistake



GENOME = ''
TMS_THRESH_GREEN = 0.7
TMS_THRESH_YELLOW = 0.5
EXCELLENT_EVAL = 1e-30
E_VAL_GREEN = 1e-10
E_VAL_YELLOW = 1e-3 
Q_COV_THRESH = 80
S_COV_THRESH = 80
LOW_COV = 50
AUTO_RED = 20
Membraneprotein_threshold=3
Hit_TMS_Diff= 2
MIN_PROTEINS = 0.5
MIN_HIT_THRESH = 100
FEATURE_TABLE = ''
GBFF_FILE = ''

superfamilies = ['1.A.17', '2.A.7', '1.E.11', '1.E.36', '1.C.12', '1.A.72', '2.A.29', '2.A.66','3.A.3', '1.C.41', '8.A.180', '2.A.6']
tcdbSystems = {}
red = {}
green = {}
yellow = {}
hmmtop_df = pd.DataFrame(columns=['Hit_tcid', 'Hit_xid', 'Hit_n_TMS','Match_length'])
geneFusions={}
df = None
id_dict = {}
tcid_assignments = {}
columns = ['Hit_tcid', 'Hit_xid', '#Query_id','Match_length','e-value','%_identity','Query_Length','Hit_Length','Q_start',
      'Q_end','S_start','S_end','Query_Coverage','Hit_Coverage','Query_n_TMS','Hit_n_TMS','TM_Overlap_Score','Family_Abrv'
      ,'Predicted_Substrate','Query_Pfam','Subject_Pfam', 'isFusion', 'Missing_components', 'Initial_decision']
green_df = pd.DataFrame(columns=columns)
yellow_df = pd.DataFrame(columns=columns)
red_df = pd.DataFrame(columns=columns)
family_assignments = {}
proteins_found = {}
tcdb_pfams = {}
genome_pfams = {}
qids = []
multi_system_assignments = {}
missing = []
missing_info = {}
queries = {}
fam_not_found = []
move_to_green = []
move_to_red = []
completed_info = {}
written = []
query_assignments = {}
tcid_queries = {}
matches = {}
new_genomic_context_matches = {}
multicomp_queries = []
curr_hit = {}
mistakes = {}
missing_overlap_dict = {}
missing_tms_dict = {}
already_swapped = {}
new_matches = {}
operons = {}
tcids_genome_change = []
mcs = {}
# assign genome pfams in function that ids missing components isMultiComp(row)
def assign_tcdb_pfams():
    # this file never changes so if it exists we would just use that
    if not os.path.exists(GENOME + '/analysis/tcdb_pfam.out'):
        cmd = f"bzcat $TCDBPFAM | grep -v '#' | perl -pe 's/ +/\t/g;' | cut -f 2,4 > {GENOME}/analysis/tcdb_pfam.out"
        os.system(cmd)

    with open(GENOME + '/analysis/tcdb_pfam.out', 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        content = line.split('\t')
        if content[1] in tcdb_pfams:
            tcdb_pfams[content[1]].append(content[0].split('.')[0])
        else:
            tcdb_pfams[content[1]] = [content[0].split('.')[0]]


# This function gives each family an assignment to be used when identifying if a family is in green
# creates the family_assignents dictionary that contains all info for missing component analysis
def assign_family(row, color):
    row_index = 0
    row_as_list = row.iloc[row_index].to_dict()
    tcid = row['Hit_tcid'].tolist()[0]
    xid = row_as_list['Hit_xid']
    # if the tcid is not found, we want to instantiate the array
    if tcid not in proteins_found:
        proteins_found[tcid] = []
    # if the protein is already in use in a family, we dont need to assign it again
    if xid in proteins_found[tcid]:
        return
    else:
        proteins_found[tcid].append(xid)
    temp = tcid.split('.')
    # define a family either on the first 3 characters of tcid or first 4
    fam = temp[0] + '.' + temp[1] + '.' + temp[2]
    if fam in superfamilies:
        fam = temp[0] + '.' + temp[1] + '.' + temp[2] + '.' + temp[3]
    # assign families to colors
    if color not in family_assignments:
        family_assignments[color] = {}
        if fam not in family_assignments[color]:
            family_assignments[color][fam] = {}
        if tcid in family_assignments[color][fam]:
            family_assignments[color][fam][tcid].append(row_as_list)
        else:
            family_assignments[color][fam][tcid] = [row_as_list]
    
    else:
        if fam not in family_assignments[color]:
            family_assignments[color][fam] = {}
        if tcid in family_assignments[color][fam]:
            family_assignments[color][fam][tcid].append(row_as_list)
        else:
            family_assignments[color][fam][tcid] = [row_as_list]

'''
def create_id_dict(df):
    if df == None:
        return {}
    id_df = df[['Hit_tcid', 'Hit_xid', '#Query_id','Match_length','e-value','%_identity','Query_Length','Hit_Length','Q_start',
     'Q_end','S_start','S_end','Query_Coverage','Hit_Coverage','Query_n_TMS','Hit_n_TMS','TM_Overlap_Score','Family_Abrv'
     ,'Predicted_Substrate','Query_Pfam','Subject_Pfam']]
    
    id_dict = {}
    for index, row in id_df.iterrows():
        if row['Hit_tcid'] in id_dict:
            id_dict[row['Hit_tcid']].append(row)
        else:
            id_dict[row['Hit_tcid']] = [row]
    return id_dict
'''

# there is an error in gblast's calculation of tms overlap so this uses parseXML to calculate correct values to be changed in df
def adjustOverlapScore(df):
    print(GENOME)
    # this calls the parseXML script to adjust overlap scores of tmss
    overlap_dict = parse(GENOME + '/hmmtop.db', GENOME + '/xml/' ,GENOME + '/results.tsv', 8)
    score_dict = {}
    # corrects the overlap score in the table itself
    for k in overlap_dict:
        score_dict[k] = overlap_dict[k]['alignedTMS']
    for index, row in df.iterrows():
        if row['#Query_id'] in score_dict:
            df.at[index, 'TM_Overlap_Score'] = score_dict[row['#Query_id']]
        else:
            df.at[index, 'TM_Overlap_Score'] = 0
    print('finished')

#This function will fill the dictionary tcdbSystems with data.
def parseTCDBcontent(input):
    for line in input:
        if(">" in line):
            first = line.split("-")[0][1:]
            if(first in tcdbSystems):
                tcdbSystems.get(first).append(line.strip(">\n"))
            else:
                tcdbSystems[first] = [line.strip(">\n")]
#This is a helper function for finding common pFam domains and can be used to check if a value is a float
def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

# this function returns the query coverage of the given row
def qCoverage(row):
    return float(row.get(key='Query_Coverage'))

# this function returns the target coverage of the given row
def hCoverage(row):
    return float(row.get(key='Hit_Coverage'))

# this function returns the evalue of the given row
def eVal(row):
    return float(row.get(key='e-value'))

# this function returns the common pfam domains of the given row
def pfamDoms(row):
    doms = []
    if isfloat(row.get(key='Query_Pfam')):
        return doms
    elif isfloat(row.get(key='Subject_Pfam')):
        return doms

    q_pfam = row.get(key='Query_Pfam').split(',')
    s_pfam = row.get(key='Subject_Pfam').split(',')
    
    # if there are no pfams on either protein mark as NA
    if len(q_pfam) == 0 and len(s_pfam) == 0:
        return ['N/A']
    for q_domain in q_pfam:
        for s_domain in s_pfam:
            if q_domain == s_domain:
                doms.append(q_domain)
    common_doms = [*set(doms)]
    return common_doms

# this function makes a decision based on the tms overlap of the proteins
def tm_overlap_decision(row):
    q_tmss = int(row.get(key='Query_n_TMS'))
    t_tmss = int(row.get(key='Hit_n_TMS'))
    # if there are no tms regions in the tcdb protein tms regions are irrelevant in decision making
    if t_tmss == 0 and q_tmss == 0:
        return True
    # calculates the percent overlap of tms regions relative to tcdb proteins
    def tmoverlap_percent(row):
        overlap_percent = int(row.get(key='TM_Overlap_Score')) / t_tmss
        return overlap_percent

    # if there are under 3 tmss then theres a possibility of there being mishits
    if t_tmss > 3:
        # if there is great overlap return true
        if tmoverlap_percent(row) >= TMS_THRESH_GREEN:
            return True
    else:
        # could potentially have a fusion candid so should be placed into yellow
        if q_tmss >= 1 and t_tmss >= 1:
            return True
    
    return False
        
# This function calculates the tms overlap percentage between the 2 proteins relative to the tcdb protein
def overlap_percent(row):
    q_tmss = int(row.get(key='Query_n_TMS'))
    t_tmss = int(row.get(key='Hit_n_TMS'))
    if q_tmss == 0 or t_tmss == 0:
        return 0
    
    overlap_percent = int(row.get(key='TM_Overlap_Score')) / max(q_tmss, t_tmss)
    return overlap_percent

# this makes a decision based on the coverage of the fused proteins together
def assign_if_fusion(fusions, tcdb_tms):
    # need to check that the combined tms count from all fusions meets threshold if ritvik can add to output
    tot_cov = 0
    # print(fusions)
    if len(fusions) == 0:
        return 'red'
    tot_cov = int(fusions[0]['send']) - int(fusions[0]['sstart'])
    for i in range(1, len(fusions)):
        # caluculate overlap
        overlap = 0
        if fusions[i]['sstart'] < fusions[i-1]['send']:
            overlap = fusions[i-1]['send'] - fusions[i]['sstart']
        tot_cov += fusions[i]['send'] - fusions[i]['sstart'] - overlap
        if fusions[i]['send'] < fusions[i-1]['send']:
            return 'red'
        if fusions[i]['sstart'] == fusions[i-1]['sstart']:
            return 'red'

    # if there are no tmss assume theres 100% overlap
    if tcdb_tms == 0:
        overlap_percent_fus = 1
    else:
        net_overlap = fusion_TMS_count(fusions)
        overlap_percent_fus = net_overlap / tcdb_tms
        

    tot_cov = float(tot_cov / int(fusions[0]['hit_length'])) * 100
    
    # if the overall coverage covers most of the protein and there is good tms overlap its a green
    if tot_cov >= S_COV_THRESH and overlap_percent_fus >= TMS_THRESH_GREEN:
        return 'green'
    elif tot_cov >= LOW_COV and overlap_percent_fus >= TMS_THRESH_YELLOW:
        return 'yellow'
    else:
        return 'red'
    
# this is the decision making algorithm based on all essential values from blast results
def make_decision(row, fusions, tcdb_tms):
    eval = eVal(row)
    qcov = qCoverage(row)
    scov = hCoverage(row)
    pfam_doms = pfamDoms(row)
    tmoverlap = overlap_percent(row)
    # automatic green if e-val is 0 or extremely good with okay coverage and has some tms coverage
    if (eval == 0.0 or eval <= EXCELLENT_EVAL) and (qcov > LOW_COV and scov > LOW_COV) and tcdb_tms > 0:
        return 'green'
    elif (eval == 0.0 or eval <= EXCELLENT_EVAL) and (qcov > Q_COV_THRESH and scov > S_COV_THRESH):
        return 'green'
    
    # if coverage is less than a very low number it is impossible for it to be a good hit
    if qcov <= AUTO_RED or scov <= AUTO_RED and len(fusions) == 0:
        return 'red'
    # if there are less than 3 tms, there is a possiblity of mischaracterizing tms regions
    if tcdb_tms > 3:
        # great e val, there are common pfam domains, good cov and high tms overlap means good hit
        if eval <= E_VAL_GREEN and len(pfam_doms) > 0 and qcov >= Q_COV_THRESH and scov >= S_COV_THRESH and tmoverlap >= TMS_THRESH_GREEN:
            return 'green'
        # great e val, common doms, has good coverage or there is a high tms overlap means good hit
        if eval <= E_VAL_GREEN and len(pfam_doms) > 0 and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) or tmoverlap >= TMS_THRESH_GREEN):
            return 'green'
        # great e value, has good coverage or there is a high tms overlap means good hit
        if eval <= E_VAL_GREEN and len(pfam_doms) > 0 and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) or tmoverlap >= TMS_THRESH_GREEN):
            return 'green'
        # okay e value, has good coverage and there is a high tms overlap means good hit
        if eval <= E_VAL_YELLOW and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) and tm_overlap_decision(row) == True):
            return 'green'
        # great e-val no matching pfam domains but has goood coverage or good tmoverlap its a good hit
        if eval <= E_VAL_GREEN and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) or tmoverlap >= TMS_THRESH_GREEN):
            return 'green'
        # okay eval, full tms overlap w tcdb, matching pfam doms, even if coverage isnt great is a green
        if eval <= E_VAL_GREEN and ((qcov >= LOW_COV and scov >= LOW_COV) and tm_overlap_decision(row) == True) and len(pfam_doms) > 0:
            return 'green'
        # okay eval, good tms overlap and ok coverage is a green
        if eval <= E_VAL_YELLOW and (qcov >= LOW_COV and scov >= LOW_COV) and tm_overlap_decision(row) == True:
            return 'green'
        # ok eval with some tms overlap should be green
        if eval <= E_VAL_YELLOW and (qcov >= AUTO_RED and scov >= AUTO_RED) and tmoverlap >= TMS_THRESH_YELLOW:
            return 'green'
        # if there is a fusion found we need to do further analysis
        if len(fusions) > 0:
            return 'fusion found'
        
        # ok e-val, no common doms, has good coverage or there is a good tms overlap means yellow hit
        if eval <= E_VAL_YELLOW and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) or tmoverlap >= TMS_THRESH_YELLOW):
            return 'yellow'
        # okay e val, and com dom tms overlap is good and the coverage is not too terrible > 50% its a yellow hit
        if eval <= E_VAL_YELLOW and len(pfam_doms) > 0 and (qcov > LOW_COV or scov > LOW_COV) and tmoverlap >= TMS_THRESH_YELLOW:
            return 'yellow'
        # okay e val, common domains, with good tms overlap even if coverage is low is yellow
        if eval <= E_VAL_YELLOW and len(pfam_doms) > 0 and tm_overlap_decision(row) == True:
            return 'yellow'
        # okay eval, good tms overlap is yellow could have fusion not listed
        if eval <= E_VAL_YELLOW and ((qcov >= LOW_COV and scov >= LOW_COV) and tm_overlap_decision(row) == True):
            return 'yellow'
        # if low e-val, has low coverage and there are low tms overlap
        if eval <= E_VAL_YELLOW and (qcov < LOW_COV or scov < LOW_COV) and tmoverlap < TMS_THRESH_YELLOW:
            return 'red'
    else: # considering possibility of mischaracterization of tms
        # is a fusion candidate, has great e value and has good cov is good hit
        if eval <= E_VAL_GREEN and (qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) and len(pfam_doms) > 0:
            return 'green'
        # good e value, common doms match, good cov is good hit
        if eval <= E_VAL_GREEN and len(pfam_doms) > 0 and (qcov >= Q_COV_THRESH and scov >= S_COV_THRESH):
            return 'green'
        # great e-val no matching pfam domains but has goood coverage or good tmoverlap its a good hit
        if eval <= E_VAL_YELLOW and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) and tm_overlap_decision(row) == True):
            return 'green'
        # ok evalue with matching pfam domains should be green
        if eval <= E_VAL_YELLOW and (qcov >= LOW_COV and scov >= LOW_COV) and len(pfam_doms) > 0:
            return 'green'
        # if there is a fusion found we need to do further analysis
        if len(fusions) > 0:
            return 'fusion found'
        
        # has ok e value, common doms match, good cov is good hit
        if eval <= E_VAL_YELLOW and len(pfam_doms) > 0 and (qcov >= LOW_COV and scov >= LOW_COV):
            return 'yellow'
        # is fusion, ok e value and good coverage but not good tms overlap
        if eval <= E_VAL_YELLOW and qcov >= Q_COV_THRESH and scov >= S_COV_THRESH:
            return 'yellow'
        # great e val, and com dom but not good coverage is yellow hit
        if eval <= E_VAL_GREEN and len(pfam_doms) > 0:
            return 'yellow'
        # okay e val, and com dom but not good coverage is yellow hit
        if eval <= E_VAL_YELLOW and len(pfam_doms) > 0:
            return 'yellow'
        # okay e val, no common doms, good tms overlap
        if eval <= E_VAL_YELLOW and tmoverlap >= TMS_THRESH_GREEN:
            return 'yellow'
        
    return 'red'

# based on all factors and using make_decision as a helper, it makes a final decision on the hit
def final_decision(row, fusions, tcdb_tms):
    eval = eVal(row)
    qcov = qCoverage(row)
    scov = hCoverage(row)
    pfam_doms = pfamDoms(row)
    tmoverlap = overlap_percent(row)
    decision = make_decision(row, fusions, tcdb_tms)
    # if theres a fusion found we have to decide differently based on overall coverage
    if decision == 'fusion found':
        # if assign_if_fusion returns red its a false fusion or unlikely
        if assign_if_fusion(fusions, tcdb_tms) == 'red':
            return make_decision(eval, qcov, scov, pfam_doms, [], tmoverlap, tcdb_tms)
        else:
            return assign_if_fusion(fusions, tcdb_tms)
    else:
        return decision


# This function will check the tcdbSystems dictionary for the tcid given 
def isSingleComp(row):
    tcid = row["Hit_tcid"]
    tc_arr = tcdbSystems.get(tcid)
    if(len(tc_arr) == 1):
        return True
    else:
        return False
    
def categorizeSingleComp(row):
    row_to_write = row.tolist()
    sortedGeneArr = []
    fusions = ''
    tcid = row['Hit_tcid']
    tcdb_tms = row['Hit_n_TMS']    
    fus_list = []
    fusion_results= geneFusions[row["Hit_tcid"] + "-" + row["Hit_xid"]]
    
    # check to see if potential fusion
    if(len(fusion_results) !=1):
        sortedGeneArr = sorted(fusion_results, key=lambda x: x['sstart'])
        if(len(isFusion(sortedGeneArr))!=0):
            sortedGeneArr = isFusion(sortedGeneArr)
            # for any potential fusion append to fusions list
            for k in range(0, len(sortedGeneArr)):
                fusions += sortedGeneArr[k]['query'] + ' '
            row_to_write.append(fusions)
        else:
            sortedGeneArr = []
            row_to_write.append(fusions)
    else:
        row_to_write.append(fusions)
    if assign_if_fusion(sortedGeneArr, tcdb_tms) == 'red':
        sortedGeneArr = []
    
    if len(fusions) > 0:
        fus_list = fusions.split(' ')
        if row['#Query_id'] not in fus_list:
            row_to_write.pop()
            row_to_write.append('')
    fus_color = ''

    if len(fusions) != 0:
        fus_color = assign_if_fusion(sortedGeneArr, tcdb_tms)

    # if there is a copy we need to assign it to the same color as the other copies
    if tcid in tcid_assignments:
        return (tcid_assignments[tcid], row_to_write)
    
    decision = final_decision(row, sortedGeneArr, tcdb_tms)
    # based on the decision send it to the appropriate file
    if decision == 'green':
        green[row.get(key='#Query_id')] = row_to_write
        tcid_assignments[row.get(key='Hit_tcid')] = 'Green'
        queries[row['#Query_id']] = 'Green'
        return ('Green', row_to_write)
    elif decision == 'yellow':
        yellow[row.get(key='#Query_id')] = row_to_write
        tcid_assignments[row.get(key='Hit_tcid')] = 'Yellow'
        queries[row['#Query_id']] = 'Yellow'
        return ('Yellow', row_to_write)
    elif decision == 'red':
        red[row.get(key='#Query_id')] = row_to_write
        tcid_assignments[row.get(key='Hit_tcid')] = 'Red'
        queries[row['#Query_id']] = 'Red'
        return ('Red', row_to_write)

    

# This command executes bash commands on the command line
def execute_command(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    process.wait()
    if process.returncode == 0:
        print(f"command'{command}'success")
    else:
        print(f"command'{command}'fail")

# returns empty list when there are no membrane proteins found. 
# otherwise returns list of membrane protein accessions.
def FindMembraneProtein(row,df):
    
    if isSingleComp(row): ##get rid of single comp system first
        return ([])
    tcid = row["Hit_tcid"]
    tc_arr = tcdbSystems.get(tcid)##This contains all the proteins with their system name

    tc_all_arr= [arr.split("-")[1] for arr in tc_arr] ##All the proteins that TCDB provide
    
    ##if tcid.split(".")[1]=="B":##and 
    if tcid[0:3]=="1.B":##If there is a beta barrel, we assume they are membrane proteins
        return (tc_all_arr)
    tc_filter_arr= list(filter(lambda x: (len(df.query(f"Hit_xid=='{x}' and Hit_tcid=='{tcid}'"))) !=0, tc_all_arr))
    ##tc_filter_arr includes the proteins that showed up in the gblast result
    tc_Notfind_arr=list(set(tc_all_arr)-set(tc_filter_arr))##This is the missing proteins in actual result
    find_df = pd.concat([df.query(f"Hit_xid=='{arr}' and Hit_tcid=='{tcid}'") for arr in tc_all_arr])
    ##This is the whole line for tc_filter_arr in gblast result
    unfind_df = pd.concat([hmmtop_df.query(f"Hit_xid=='{arr}' and Hit_tcid=='{tcid}'") for arr in tc_all_arr])
    ##This is the whole line for tc_all_arr in hmmtop file
    if find_df['Hit_n_TMS'].nunique() == 1 and unfind_df['Hit_n_TMS'].nunique() == 1 :
        ##If all the proteins have same Hit TMS in both files, we say they are membrane proteins
        if str(find_df['Hit_n_TMS'].unique()[0]) == str(unfind_df['Hit_n_TMS'].unique()[0]):
           
            return (tc_all_arr)
    Found_arr_gblast = find_df[(find_df['Hit_n_TMS'] >= Membraneprotein_threshold) & (abs(find_df['Hit_n_TMS'] - find_df['Query_n_TMS']) <= Hit_TMS_Diff)]["Hit_xid"].tolist()
    ##If there are proteins that has Hit TMS >=3 and difference with its query TMS <=2, we say it's membrane proteins
    tcdb_arr= Found_arr_gblast +tc_Notfind_arr##This arr contains the possible output proteins 
    if len(tcdb_arr)==0:
        return([])
    tcdb_df = pd.concat([hmmtop_df.query(f"Hit_xid=='{arr}' and Hit_tcid=='{tcid}'") for arr in tcdb_arr])
    Final_return = tcdb_df[(tcdb_df['Hit_n_TMS']>= Membraneprotein_threshold)]  
    ##print(tc_all_arr)
    ##print(find_df)
    ##print(unfind_df)
    ##print(Final_return["Hit_xid"].tolist())
    # return((Final_return['Hit_tcid'] + '-' + Final_return["Hit_xid"]).tolist())  

    return Final_return.to_dict()
    ##final_df = pd.concat([df.query(f"Hit_xid=='{arr}'") for arr in tc_filter_arr])

# determines if membrane protein dict is empty
def isEmpty(dict):
    if type(dict) is list:
        return len(dict) == 0
    count = 0
    for key in dict:
        if not bool(dict[key]):
            count += 1
    
    return count == 4

# if there are missing components we need to to take a different decision determining the importance of the components that are missing
def multicomp_decision(tcdb_proteins, tcid, mem_dict, tc_filter_arr, missing_comps, has_file):
    
    #find_pfams(missing_comps, tcid)
    protein_type = tcid.split('.')[0] + '.' + tcid.split('.')[1]
    max_tms = 0
   
    tcid = tcdb_proteins[0].split('-')[0]
    # if there are no membrane proteins and we hit the min protein threshold
    if (isEmpty(mem_dict) or len(mem_dict) == 0) and (len(tc_filter_arr)/len(tcdb_proteins)) < MIN_PROTEINS:
        return 'red'
    # if there are no membrane proteins keep in yellow for now
    elif isEmpty(mem_dict) or len(mem_dict) == 0:
        return 'yellow'
    
    membrane_protein_tms = {}
    
    # put membrane proteins into a dict sorted by number of tmss
    membrane_proteins = []
    if type(mem_dict) is list:
        membrane_proteins = mem_dict
    else:
        membrane_proteins_accessions = list(mem_dict["Hit_xid"].values())
        membrane_proteins_tcids = list(mem_dict["Hit_tcid"].values())
        membrane_tmss = list(mem_dict["Hit_n_TMS"].values())
        
        membrane_proteins = membrane_proteins_accessions
        for i in range(0, len(membrane_proteins_accessions)):
            protein = membrane_proteins_tcids[i] + '-' + membrane_proteins_accessions[i]
            membrane_protein_tms[protein] = membrane_tmss[i]
            membrane_tmss = list(mem_dict['Hit_n_TMS'].values()) 
        # sorts the dictionary based on ascending tms values
        membrane_protein_tms = {k: membrane_protein_tms[k] for k in sorted(membrane_protein_tms, key=membrane_protein_tms.get, reverse=True)}
        for key in membrane_protein_tms:
            max_tms = membrane_protein_tms[key]
            break

    # if beta barrell and most proteins are in the system should be yellow->green because ALL are membrane proteins
    if protein_type == '1.B' and abs(len(tcdb_proteins) - len(tc_filter_arr)) <= 2:
        return 'Green'
    all_mems_found = True

    # determine if all membrane proteins are found if so its a green
    if type(mem_dict) is list:
        for protein in membrane_proteins:
            if protein not in tc_filter_arr:
                all_mems_found = False
                break
    else:
        for protein in membrane_proteins:
            accession = protein
            if accession not in tc_filter_arr:
                all_mems_found = False
                break
    # if we have all the membrane proteins we can confidently say the system exists or performs the same function
    if all_mems_found == True:
        return 'Green'

    # if we dont have enough proteins under 50% usually
    if len(tcdb_proteins) > 3 and (len(tc_filter_arr)/len(tcdb_proteins)) < MIN_PROTEINS:
        return 'red'
    
    genome_membrane_accessions = tc_filter_arr

    num_membrane_proteins = len(genome_membrane_accessions)            
    
    # it is an instant red if NONE of the membrane proteins are found
    if num_membrane_proteins == 0:
        return 'red'
    
    # check to see if the tms with equal or most tms is found
    for protein in genome_membrane_accessions:
        key = tcid + '-' + protein
        
        # if we find the most significant membrane protein it is worth considering for green
        if key in membrane_protein_tms:
            if max_tms - membrane_protein_tms[key] < 2:
                return 'yellow'
        
    return 'red'


# Goal of this function is to run pfam.sh on ALL missing components for the second iteration
def run_query_pfam(has_file):
    genome = ''
    # extract genome name to be used to find blast results
    for comp in GENOME.split('/'):
        if 'GCF' in comp:
            genome = comp
    if not has_file:
        # write query ids to an output file to be used later
       with open(GENOME + '/analysis/query_ids.txt', 'w') as f:
            
            for q in missing_info:
                if len(missing_info[q]) == 0:
                    continue
                qid = ''
                for i in missing_info[q]['query_id'].split('/'):
                    if '.xml' in i:
                        qid = i.split('.')[0] + '.' + i.split('.')[1]
                f.write(qid + '\n') 

    if not os.path.exists(GENOME + '/analysis/query_pfam.out'):
        if not os.path.isdir(GENOME + '/analysis/blastdb'):
            create_blastdb = f'gunzip -c /ResearchData/Microbiome/Assemblies/{genome}/*.faa.gz | makeblastdb -dbtype prot -input_type fasta -parse_seqids  -hash_index -title {genome} -out blastdb/{genome} -in -; mkdir {GENOME}/analysis/blastdb; mv blastdb/* {GENOME}/analysis/blastdb'
            print(create_blastdb)
            os.system(create_blastdb)

        cmd1 = f'blastdbcmd -db {GENOME}/analysis/blastdb/{genome} -entry_batch {GENOME}/analysis/query_ids.txt -target_only > {GENOME}/analysis/queries.faa'
        print(cmd1)
        os.system(cmd1)
        # get command to run pfam with a list of query sequences
        cmd2 = f'run_pfam.sh {GENOME}/analysis/queries.faa {GENOME}/analysis/query_pfam.out'
        print(cmd2)
        os.system(cmd2)
        cmd3 = f"cat {GENOME}/analysis/query_pfam.out | grep -v '#' | perl -pe 's/ +/\\t/g;' | cut -f 2,4 > {GENOME}/analysis/query_pfam_filtered.out"
        print(cmd3)
        os.system(cmd3)

    with open(GENOME + '/analysis/query_pfam_filtered.out', 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        content = line.split('\t')
        if content[1] in genome_pfams:
            genome_pfams[content[1]].append(content[0].split('.')[0])
        else:
            genome_pfams[content[1]] = [content[0].split('.')[0]]

# On the first iteration we want a relative decision for multicomp systems before checking missing_comp data
def assign_multicomp_sys(row,df,input,has_file):
    tcid = row["Hit_tcid"]
    Fusion_Add=[]
    accession = row["Hit_tcid"] + "-" + row["Hit_xid"]
    Fusion_results= geneFusions[row["Hit_tcid"] + "-" + row["Hit_xid"]]
    sortedGeneArr = []
    if(len(Fusion_results) !=1):
        sortedGeneArr = sorted(Fusion_results, key=lambda x: x['sstart'])
        if len(isFusion(sortedGeneArr)) != 0:
            Fusion_Add=[x["query"] for x in isFusion(sortedGeneArr)]
        else:
            sortedGeneArr = []
    tc_arr = tcdbSystems.get(tcid)
    tc_all_arr= [arr.split("-")[1] for arr in tc_arr] ##All the proteins that TCDB provide
    tc_filter_arr= list(filter(lambda x: (len(df.query(f"Hit_xid=='{x}' and Hit_tcid=='{tcid}'"))) !=0, tc_all_arr))
    tc_missing_arr= list(set(tc_all_arr) - set(tc_filter_arr))
    tc_id_acc = row["Hit_tcid"] + "-" + row["Hit_xid"]
    
    query = row['#Query_id']
    mcs[tc_id_acc] = query
    if tcid in tcid_queries:
        tcid_queries[tcid].append(query)
    else:
        tcid_queries[tcid] = [query]
    matches[query] = row["Hit_tcid"] + "-" + row["Hit_xid"] 
    
# This function gets the multicomp decision in the first iteration
def isMultiComp(row,df,input,has_file):
    if isSingleComp(row):
        return (False)
    ##get rid of single comp system first
    
    tcid = row["Hit_tcid"]
    Fusion_Add=[]
    accession = row["Hit_tcid"] + "-" + row["Hit_xid"]
    Fusion_results= geneFusions[row["Hit_tcid"] + "-" + row["Hit_xid"]]
    sortedGeneArr = []
    if(len(Fusion_results) !=1):
        sortedGeneArr = sorted(Fusion_results, key=lambda x: x['sstart'])
        if len(isFusion(sortedGeneArr)) != 0:
            Fusion_Add=[x["query"] for x in isFusion(sortedGeneArr)]
        else:
            sortedGeneArr = []
    tc_arr = tcdbSystems.get(tcid)
    tc_all_arr= [arr.split("-")[1] for arr in tc_arr] ##All the proteins that TCDB provide
    tc_filter_arr= list(filter(lambda x: (len(df.query(f"Hit_xid=='{x}' and Hit_tcid=='{tcid}'"))) !=0, tc_all_arr))
    tc_missing_arr= list(set(tc_all_arr) - set(tc_filter_arr))
    tc_id_acc = row["Hit_tcid"] + "-" + row["Hit_xid"]
    
    fusions = '' 
    qid = row['#Query_id']
    
    if len(Fusion_Add) > 0:
        if qid not in Fusion_Add:
            Fusion_Add = []       
 
    protein_type = row['Hit_tcid'].split('.')[0] + '.' + row['Hit_tcid'].split('.')[1]
    tcdb_tms = row['Hit_n_TMS']
    

    # if the tcid has already been assigned, assign this protein to the same file
    if row['Hit_tcid'] in tcid_assignments:
        if tcid_assignments[row['Hit_tcid']] == 'Green':
            matches[row['#Query_id']] = row['Hit_tcid'] + '-' + row['Hit_xid']
            curr_hit[tc_id_acc] = row['#Query_id']
            if len(tc_missing_arr) > 0:
                for missing in tc_missing_arr:
                    if len(missing) > 0: 
                        missing_info[row['Hit_tcid'] + '-' + missing] = {}
 
        queries[row['#Query_id']] = tcid_assignments[row['Hit_tcid']]
        tcid_queries[row['Hit_tcid']].append(row['#Query_id'])
        return({"color": tcid_assignments[row['Hit_tcid']],
                "Found_proteins":tc_filter_arr,
                "All_proteins":tc_arr,
                'Missing_proteins':tc_missing_arr,
                "Fusion_results":Fusion_Add,
                "isFusion":len(Fusion_Add)>0,
                "Initial_decision": tcid_assignments[row['Hit_tcid']]})

# This redecides mcs after looking at genomic context
def isMultiCompSecondIteration(row,df,input,has_file):
    if isSingleComp(row):
        return (False)
    ##get rid of single comp system first
    
    tcid = row["Hit_tcid"]
    Fusion_Add=[]
    accession = row["Hit_tcid"] + "-" + row["Hit_xid"]
    Fusion_results= geneFusions[row["Hit_tcid"] + "-" + row["Hit_xid"]]
    sortedGeneArr = []
    if(len(Fusion_results) !=1):
        sortedGeneArr = sorted(Fusion_results, key=lambda x: x['sstart'])
        if len(isFusion(sortedGeneArr)) != 0:
            Fusion_Add=[x["query"] for x in isFusion(sortedGeneArr)]
        else:
            sortedGeneArr = []
    tc_arr = tcdbSystems.get(tcid)
    tc_all_arr= [arr.split("-")[1] for arr in tc_arr] ##All the proteins that TCDB provide
    tc_filter_arr= list(filter(lambda x: (len(df.query(f"Hit_xid=='{x}' and Hit_tcid=='{tcid}'"))) !=0, tc_all_arr))
    tc_missing_arr= list(set(tc_all_arr) - set(tc_filter_arr))
    tc_id_acc = row["Hit_tcid"] + "-" + row["Hit_xid"]
    
    fusions = '' 
    qid = row['#Query_id']
    
    if len(Fusion_Add) > 0:
        if qid not in Fusion_Add:
            Fusion_Add = []       

    protein_type = row['Hit_tcid'].split('.')[0] + '.' + row['Hit_tcid'].split('.')[1]
    tcdb_tms = row['Hit_n_TMS']
    
    mismatched_query = ''
    # if the system has been determined to be an operon we have to find a way to complete the system with the right proteins
    if tcid in operons and operons[tcid] == 'possible operon':
        MembraneProteins= FindMembraneProtein(row, df)
        decision = multicomp_decision(tc_arr, row['Hit_tcid'], MembraneProteins, tc_filter_arr, tc_missing_arr, has_file)

        if tcid in mistakes and decision != 'red':
            
            for m in range(0, len(mistakes[tcid])):
                
                mismatched_query = mistakes[tcid][m]
                
                
                corrected_protein_assignment = matches[mismatched_query]
                if len(corrected_protein_assignment) == 0:
                    continue
                # should be corrected_protein is tcid-component with correct query
                # that current assignment in row is the wrong query
                # we need to remove that entire row and create a new one with new values
                
                protein_to_remove = corrected_protein_assignment.split('-')[1]
                if protein_to_remove in tc_filter_arr:
                    tc_filter_arr.remove(protein_to_remove)
                tc_missing_arr.append(protein_to_remove)
                new_matches[corrected_protein_assignment] = mismatched_query
                # add queries and tcid to dict and figure out if we can just automatically classify in green if part of operon according to genomic context
    # If the system has been previously assigned to a file just assign the protein to the same file
    if row['Hit_tcid'] in tcid_assignments:
        if tcid_assignments[row['Hit_tcid']] == 'Green':
            matches[row['#Query_id']] = row['Hit_tcid'] + '-' + row['Hit_xid']
            curr_hit[tc_id_acc] = row['#Query_id']
            if len(tc_missing_arr) > 0:
                for missing in tc_missing_arr:
                    if len(missing) > 0: 
                        missing_info[row['Hit_tcid'] + '-' + missing] = {}

        queries[row['#Query_id']] = tcid_assignments[row['Hit_tcid']]
        tcid_queries[row['Hit_tcid']].append(row['#Query_id'])
        return({"color": tcid_assignments[row['Hit_tcid']],
                "Found_proteins":tc_filter_arr,
                "All_proteins":tc_arr,
                'Missing_proteins':tc_missing_arr,
                "Fusion_results":Fusion_Add,
                "isFusion":len(Fusion_Add)>0,
                "Initial_decision": tcid_assignments[row['Hit_tcid']]})

    
    if(set(tc_all_arr)==set(tc_filter_arr)):##If all the proteins in that system can be found, then green
        #print(tc_arr)
        #print("green",tc_filter_arr,"Fusion Results:",Fusion_Add)
        #if tcid == '3.A.1.12.4':
        #    print(eVal(row) <= E_VAL_GREEN and qCoverage(row) >= Q_COV_THRESH and hCoverage(row) >= S_COV_THRESH and len(pfamDoms(row)) != 0 and tm_overlap_decision(row) == True)
        #if(eVal(row) <= E_VAL_GREEN and qCoverage(row) >= Q_COV_THRESH and hCoverage(row) >= S_COV_THRESH and len(pfamDoms(row)) != 0 and tm_overlap_decision(row) == True):
        if(eVal(row) <= E_VAL_GREEN and qCoverage(row) >= Q_COV_THRESH and hCoverage(row) >= S_COV_THRESH and len(pfamDoms(row)) != 0 and tm_overlap_decision(row) == True):    
            tcid_assignments[row['Hit_tcid']] = 'Green'
            queries[row['#Query_id']] = 'Green'
            curr_hit[tc_id_acc] = row['#Query_id']
            for protein in tc_filter_arr:
                if protein in multi_system_assignments:
                    multi_system_assignments[protein].append(tcid, 'Green')
                else:
                    multi_system_assignments[protein] = (tcid, 'Green')
            tcid_queries[row['Hit_tcid']] = [row['#Query_id']]
            matches[row['#Query_id']] = row['Hit_tcid'] + '-' + row['Hit_xid']
            return({"color":"Green",
                    "Found_proteins":tc_filter_arr,
                    "All_proteins":tc_arr,
                    'Missing_proteins':tc_missing_arr,
                    "Fusion_results":Fusion_Add,
                    "isFusion":len(Fusion_Add)>0,
                    "Initial_decision": 'green'})
    
    # if this is a new match determined by genomic context has to become green
    if tc_id_acc in new_matches:
        matches[new_matches[tc_id_acc]] = row['Hit_tcid'] + '-' + row['Hit_xid']
        curr_hit[tc_id_acc] = new_matches[tc_id_acc]
        if len(tc_missing_arr) > 0:
            for missing in tc_missing_arr:
                if len(missing) > 0:
                    missing_info[row['Hit_tcid'] + '-' + missing] = {}
        
        tcid_assignments[row['Hit_tcid']] = 'Green'
        queries[new_matches[tc_id_acc]] = tcid_assignments[row['Hit_tcid']]
        tcid_queries[row['Hit_tcid']] = [new_matches[tc_id_acc]]
        return({"color":"Green",
                    "Found_proteins":tc_filter_arr,
                    "All_proteins":tc_arr,
                    'Missing_proteins':tc_missing_arr,
                    "Fusion_results":Fusion_Add,
                    "isFusion":len(Fusion_Add)>0,
                    "Initial_decision": 'green'})
    # get list of membrane proteins in the system
    MembraneProteins= FindMembraneProtein(row, df)
    count_mem = 0
    total_tcmem = 0
    # count the number of membrane proteins to be compared later
    if len(MembraneProteins) > 0 or not isEmpty(MembraneProteins):
        if type(MembraneProteins) is not list:
            mem_proteins = list(MembraneProteins['Hit_xid'].values())
            total_tcmem = len(mem_proteins)
        else:
            mem_proteins = MembraneProteins
            total_tcmem = len(MembraneProteins)

        for p in mem_proteins:
            if p in tc_filter_arr:
                count_mem += 1
    
    # Determine how essential the missing components are before assignment
    decision = multicomp_decision(tc_arr, row['Hit_tcid'], MembraneProteins, tc_filter_arr, tc_missing_arr, has_file)
    if len(tc_missing_arr) == 0 and final_decision(row, sortedGeneArr, tcdb_tms) == 'yellow':
        tcid_assignments[row['Hit_tcid']] = 'Yellow'
        for protein in tc_filter_arr:
            queries[row['#Query_id']] = 'Yellow'
            if protein in multi_system_assignments:
                multi_system_assignments[protein].append(tcid, 'Yellow')
            else:
                multi_system_assignments[protein] = (tcid, 'Yellow')
        tcid_queries[row['Hit_tcid']] = [row['#Query_id']]
        return({"color":"Yellow",
                "Found_proteins":tc_filter_arr,
                "All_proteins":tc_arr,
                'Missing_proteins':tc_missing_arr,
                "Fusion_results":Fusion_Add,
                "isFusion":len(Fusion_Add)>0,
                'Initial_decision': 'Yellow'})
    # if(input*len(tc_all_arr)<=len(tc_filter_arr)) and len(set(MembraneProteins) & set(tc_filter_arr))>0:
    if (total_tcmem > 0 and input <= float(count_mem / total_tcmem)) or (len(set(tc_filter_arr))>0 and isEmpty(MembraneProteins)):
        ##given some proteins can be found while containing the membrane proteins
        if final_decision(row, sortedGeneArr, tcdb_tms) == 'yellow' or final_decision(row, sortedGeneArr, tcdb_tms) == 'green':
            if decision == 'yellow':
                queries[row['#Query_id']] = 'Yellow'
                tcid_assignments[row['Hit_tcid']] = 'Yellow'
                for protein in tc_filter_arr:
                    if protein in multi_system_assignments:
                        multi_system_assignments[protein].append(tcid, 'Yellow')
                    else:
                        multi_system_assignments[protein] = (tcid, 'Yellow')
                tcid_queries[row['Hit_tcid']] = [row['#Query_id']]
                return({"color":"Yellow",
                        "Found_proteins":tc_filter_arr,
                        "All_proteins":tc_arr,
                        'Missing_proteins':tc_missing_arr,
                        "Fusion_results":Fusion_Add,
                        "isFusion":len(Fusion_Add)>0, 
                        'Initial_decision': 'Yellow'})
            elif decision == 'Green':
                queries[row['#Query_id']] = 'Green'
                tcid_assignments[tcid] = 'Green'
                curr_hit[tc_id_acc] = row['#Query_id']
                for protein in tc_filter_arr:
                    if protein in multi_system_assignments:
                        multi_system_assignments[protein].append(tcid, 'Green')
                    else:
                        multi_system_assignments[protein] = (tcid, 'Green')
                tcid_queries[row['Hit_tcid']] = [row['#Query_id']]
                matches[row['#Query_id']] = row['Hit_tcid']  + '-' + row['Hit_xid']
                return({"color":"Green",
                        "Found_proteins":tc_filter_arr,
                        "All_proteins":tc_arr,
                        'Missing_proteins':tc_missing_arr,
                        "Fusion_results":Fusion_Add,
                        "isFusion":len(Fusion_Add)>0,
                        'Initial_decision': 'Yellow'})

    tcid_assignments[row['Hit_tcid']] = 'Red'
    queries[row['#Query_id']] = 'Red'
    for protein in tc_filter_arr:
        if protein in multi_system_assignments:
            multi_system_assignments[protein].append(tcid, 'Red')
        else:
            multi_system_assignments[protein] = (tcid, 'Red')
    tcid_queries[row['Hit_tcid']] = [row['#Query_id']]
    return({"color":"Red",
            "Found_proteins":tc_filter_arr,
            "All_proteins":tc_arr,
            'Missing_proteins':tc_missing_arr,
            "Fusion_results":Fusion_Add,
            "isFusion":len(Fusion_Add)>0,
            'Initial_decision': 'Red'})
   
# This function checks if the replicon exists
def found_replicon(gene_ranking, tcid, replicon):
    qs = tcid_queries[tcid]
    
    for query in qs:
        if query not in gene_ranking[replicon]:
            return False
    return True

# checks if teh gene rank is within the window
def find_number_outside_range(ranking, range_limit):
    arr = ranking.keys()
    for i in range(len(arr) - 1):
        if abs(arr[i] - arr[i + 1]) > range_limit:
            return ranking[arr[i + 1]]
    return None  # If all numbers are within the range, return None

# This function completes operons by looking at genomic context
def getGenomicContext(min_distance, min_missing_comps):
    # Check if FEATURE_TABLE and GBFF_FILE exist before proceeding
    if not os.path.isfile(FEATURE_TABLE):
        raise FileNotFoundError(f"The file {FEATURE_TABLE} does not exist.")

    if not os.path.isfile(GBFF_FILE):
        raise FileNotFoundError(f"The file {GBFF_FILE} does not exist.")
    
    # commands to get important values from feature table and shape of the replicon
    cmd1 = f'zgrep CDS {FEATURE_TABLE} | cut -f 7-11 | sort -k1,1 -k2,2n'
    genomic_data = subprocess.run(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout.split('\n')

    cmd2 = f'gunzip -c {GBFF_FILE} | grep LOCUS'
    replicon_shapes = subprocess.run(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout
        
    
    shapes = {}
    for replicon in replicon_shapes.split('\n'):
        temp = replicon.split(' ')
        #print(temp)
        if temp[0] == '':
            break

        line = [item for item in temp if item != '']
        
        shapes[line[1]] = line[5]
        
        
    replicons = {}
    gene_maps = {}
    # replicons is a dict that shows feature table accessible by replicon and gene maps is the same but accessed by gene/protein id
    for line in genomic_data:
        if line == '':
            continue
        line_split = line.split('\t')
        gene = line_split[4]
            
        replicon = line_split[0]
        if replicon not in replicons:
            replicons[replicon] = [{'start': line_split[1], 'end': line_split[2], 'protein': line_split[4], 'strand': line_split[3]}]
        else:
            replicons[replicon].append({'start': line_split[1], 'end': line_split[2], 'protein': line_split[4], 'strand': line_split[3]})
            

        if gene != '':
            if gene not in gene_maps:
                gene_maps[gene] = [{'start': line_split[1], 'end': line_split[2], 'replicon': line_split[0], 'strand': line_split[3]}]
            else:
                gene_maps[gene].append({'start': line_split[1], 'end': line_split[2], 'replicon': line_split[0], 'strand': line_split[3]})

    gene_ranking_by_tcid = {} 
    gene_ranking_by_rank = {}

    gene_ranking = {}
    # this generaates a ranking of genes to measure relative distance on the plasmid
    for replicon in replicons:
        rep = replicon.split('.')[0]
        gene_ranking[rep] = {}
        gene_ranking_by_rank[rep] = {}
        gene_ranking_by_tcid[rep] = {}
        count = 0
        for i in range(0, len(replicons[replicon])):
            if replicons[replicon][i]['protein'] != '':
                gene_ranking[rep][replicons[replicon][i]['protein']] = count
                gene_ranking_by_rank[rep][count] = [replicons[replicon][i]['protein']]
                if replicons[replicon][i]['protein'] in matches:
                    gene_ranking_by_tcid[rep][count] = matches[replicons[replicon][i]['protein']]
                else:
                    gene_ranking_by_tcid[rep][count] = replicons[replicon][i]['protein']
                count += 1
        
    green_queries = matches.keys()

    # to check for other matches we need to get all good hits from xml files of multi comp systems
    all_matches = extract_all(GENOME, list(matches.keys()))
    is_operon = {}
    
    # we are going to iterate through each mcs in each replicon to complete
    for i in range(0, len(list(replicons.keys()))):
        replicon = list(replicons.keys())[i].split('.')[0]
        #tcids = [item.split('-')[0] for item in list(curr_hit.keys())]
        # obtain the hits in the system
        tcids = [item.split('-')[0] for item in list(mcs.keys())]
        mc_hits = []
        for t in tcids:
            for hit in tcdbSystems[t]:
                mc_hits.append(hit)
        
        # create a window around each hit to look for operons
        for hit in mc_hits:
            allg = 'not operon'
            ref_prot = hit
            top_tcid = hit.split('-')[0]
            if ref_prot not in mcs:
                continue
            ref_query = mcs[ref_prot]
            if ref_query not in gene_ranking[replicon]:
                continue
            if i >= len(gene_maps[ref_query]):
                continue
            ref_strand = gene_maps[ref_query][i]['strand']
            ref_start = gene_maps[ref_query][i]['start']
            ref_end = gene_maps[ref_query][i]['end']
            ref_rank = gene_ranking[replicon][ref_query]


            
            window = []
            if shapes[replicon] == 'linear':
                # TODO: change 5 to a threshold or something this is just for testing
                if ref_rank -  min_distance <= 0:
                    window = [0, ref_rank + int((min_distance / 2))]
                else:
                    window = [ref_rank - min_distance, ref_rank + int((min_distance / 2))]
            # circular replicons need looking around the circle rather than linear window
            else:
                if ref_rank - (min_distance / 2) <= 0:
                    window = [len(gene_ranking[replicon]) + (ref_rank - int((min_distance / 2))), ref_rank + int((min_distance / 2))]
                else:
                    window = [ref_rank - int((min_distance / 2)), ref_rank + int((min_distance / 2))]
            # TODO: create list of number of essential proteins if all are found operon exists and the member of the fam in the window shoudl be added to mistakes and system should eb completed w that protein
            queries_to_investigate = []
            full_system = len(tcdbSystems[top_tcid])
            components = len(tcid_queries[top_tcid])
            count = 0
            # if all the components are within the window, we can safely say this is a correctly matched operon
            for n in range(window[0], window[1] + 1):
                if n >= len(gene_ranking_by_rank[replicon]):
                    continue
                protein = gene_ranking_by_rank[replicon][n][0]
                if protein in tcid_queries[top_tcid]:
                    count += 1
                    queries_to_investigate.append(protein)
            
            if (components - count) >= min_missing_comps:
                allg = 'not operon'
                continue

            found_membrane = False
            if count == components:
                found_membrane = True
                if components == full_system:
                    allg = 'possible operon'
                    continue
            
            for query in queries_to_investigate:
                    

                # for each query we are going to check if there is a possibility of a missmatch if a component is past the distance threshold
                if query not in gene_ranking[replicon]:
                    allg = 'not operon'
                    break
                data = gene_maps[query][i]
                strand = data['strand']

                rank = gene_ranking[replicon][query]
                # if a member of the family is one of the results and the protein is within the neighborhood or right next to it try and complete the operon
                temp = top_tcid.split('.')
                fam = temp[0] + '.' + temp[1] + '.' + temp[2]
                if fam in superfamilies:
                    fam = temp[0] + '.' + temp[1] + '.' + temp[2] + '.' + temp[3]

                protein_to_move = ''
                target_proteins = []
                num_protein_components = []
                found_per_system = {}
                    

                for j in range(window[0], window[1] + 1):
                     
                    if j >= len(gene_ranking_by_rank[replicon]):
                        continue

                    protein = gene_ranking_by_rank[replicon][j][0]
                    
                    if protein in matches and fam in matches[protein].split('-')[0]:
                        num_protein_components.append((len(tcdbSystems[matches[protein].split('-')[0]]), matches[protein].split('-')[0]))
                        
                        # if its already a component of the system we can just ignore it
                        if top_tcid == matches[protein].split('-')[0]:
                            if top_tcid in found_per_system:
                                found_per_system[matches[protein].split('-')[0]] += 1
                            else:
                                found_per_system[matches[protein].split('-')[0]] = 1

                            continue
                        
                        else:
                            if matches[protein].split('-')[0] in found_per_system:
                                found_per_system[matches[protein].split('-')[0]] += 1
                            else:
                                found_per_system[matches[protein].split('-')[0]] = 1
                        
                        target_proteins.append(protein)
                    
                    
                
                
                sorted_comps = sorted(num_protein_components, key=lambda x:x[0])
                if len(target_proteins) == 1:
                    target_tcid = sorted_comps[0][1]
                    protein_to_move = target_proteins[0]
                    # make swap in in the matches
                    if top_tcid in mistakes:
                        mistakes[top_tcid].append(protein_to_move)
                    else:
                        mistakes[top_tcid] = [protein_to_move]

                    assignment = ''
                    for h in all_matches[protein_to_move]:
                        if all_matches[protein_to_move][h].split('-')[0] == top_tcid:
                            assignment = all_matches[protein_to_move][h]
                            break
                    #print(top_tcid)
                    #print(assignment)
                    #print(protein_to_move)
                    #exit()
                    if assignment == '':
                        continue
                    matches[protein_to_move] = assignment
                    allg = 'possible operon'
                     
                    missing_info[assignment] = {}
                    tcids_genome_change.append(assignment.split('-'))       
                    already_swapped[protein_to_move] = assignment
                    '''
                    tcdbSystems[top_tcid].remove(matches[match])
                    protein_to_move = matches[protein_to_move]
                    tcdbSystems[top_tcid].append(top_tcid + '-' + protein_to_move.split('-')[1])
                    missing_info[top_tcid + '-' + protein_to_move.split('-')[1]] = {}

                    tcdbSystems[protein_to_move.split('-')[0]].remove(protein_to_move)
                    tcdbSystems[protein_to_move.split('-')[0]].append(protein_to_move.split('-')[0] + '-' + matches[match].split('-')[1])
                    missing_info[protein_to_move.split('-')[0] + '-' + matches[match].split('-')[1]] = {}
                    allg = 'possible operon'
                    '''
                elif len(target_proteins) > 1:
                    target_tcid = sorted_comps[0][1]

                    distance_dict = []
                    for p in target_proteins:
                        curr_data = gene_maps[p][i]
                        if curr_data['strand'] != ref_strand:
                            continue

                        distance_dict.append((p, abs(int(curr_data['start']) - int(ref_end))))

                    sorted_distances = sorted(distance_dict, key=lambda x:x[1])
                    
                    protein_to_move = ''
                    assignment = ''
                    for d in range(0, len(sorted_distances)):
                        protein_to_move = sorted_distances[d][0]
                    
                        if len(sorted_distances) == 0:
                            continue
                        #protein_to_move = sorted_distances[0][0]
                        if protein_to_move in already_swapped:
                            continue
                        for h in all_matches[protein_to_move]:
                            if all_matches[protein_to_move][h].split('-')[0] == target_tcid:
                                assignment = all_matches[protein_to_move][h]
                                break
                        if top_tcid in mistakes:
                            mistakes[top_tcid].append(protein_to_move)
                        else:
                            mistakes[top_tcid] = [protein_to_move]
                    
                    if assignment == '':
                        continue
                    
                    matches[protein_to_move] = assignment
                    allg = 'possible operon'
                    missing_info[assignment] = {}
                    tcids_genome_change.append(assignment.split('-'))
                    already_swapped[protein_to_move] = assignment

            operons[ref_prot.split('-')[0]] = allg

def assign_multicomp(Output_dict,Output_df_row): 
    # if len(Output_dict['Missing_proteins']) != 0:
    #     print('here')
    fusions = ''
    for fusion in Output_dict['Fusion_results']:
        fusions += fusion + ' '
    Intermediate = Output_df_row.copy()
    Intermediate['isFusion'] = fusions
    if len(Output_dict['Missing_proteins']) == 0:
        Intermediate['Missing_components'] = 'NA'
    else:
        Intermediate['Missing_components'] = str(Output_dict['Missing_proteins'])
    Intermediate['Initial_decision'] = Output_dict['Initial_decision']
    filename=f"{Output_dict['color']}.tsv"
    filename = GENOME + '/analysis/' + filename
    filemode='a' if os.path.exists(filename) else 'w' 
    #print(Intermediate)
    assign_family(Intermediate, Output_dict['color']) 
    for hit_xid in Output_dict['Missing_proteins']:
        _Intermediate = Output_df_row.copy()
        _Intermediate = _Intermediate.applymap(lambda x: 'NA')

        _Intermediate["Hit_tcid"] = Output_df_row["Hit_tcid"]
        _Intermediate["Hit_xid"] = hit_xid
        Missing_infor = hmmtop_df.loc[hmmtop_df['Hit_xid'] == hit_xid]

        _Intermediate["Match_length"] = Missing_infor["Match_length"].iloc[0] if not Missing_infor.empty else "NA"
        _Intermediate["Hit_n_TMS"] = Missing_infor["Hit_n_TMS"].iloc[0] if not Missing_infor.empty else "NA"
        _Intermediate["Missing_components"] = str(Output_dict['Missing_proteins'])
        assign_family(_Intermediate, Output_dict['color'])

def Write_multicomp(Output_dict,Output_df_row): 
    # if len(Output_dict['Missing_proteins']) != 0:
    #     print('here')
    fusions = ''
    for fusion in Output_dict['Fusion_results']:
        fusions += fusion + ' '
    Intermediate = Output_df_row.copy()
    Intermediate['isFusion'] = fusions
    if len(Output_dict['Missing_proteins']) == 0:
        Intermediate['Missing_components'] = 'NA'
    else:
        Intermediate['Missing_components'] = str(Output_dict['Missing_proteins'])
    Intermediate['Initial_decision'] = Output_dict['Initial_decision']
    filename=f"{Output_dict['color']}.tsv"
    filename = GENOME + '/analysis/' + filename
    filemode='a' if os.path.exists(filename) else 'w' 
    #print(Intermediate)
    with open(filename, mode=filemode, encoding='utf-8') as f:
        Intermediate.to_csv(f, sep='\t', header=filemode=='w', index=False)
        assign_family(Intermediate, Output_dict['color']) 
    for hit_xid in Output_dict['Missing_proteins']:
        _Intermediate = Output_df_row.copy()
        _Intermediate = _Intermediate.applymap(lambda x: 'NA')

        _Intermediate["Hit_tcid"] = Output_df_row["Hit_tcid"]
        tcid = str(_Intermediate['Hit_tcid']).split()[1].strip()
        
        _Intermediate["Hit_xid"] = hit_xid
        if (tcid + '-' + hit_xid) not in written:
            if (tcid + '-' + hit_xid) in completed_info:
                hit_info = completed_info[tcid + '-' + hit_xid]
                
                qid = ''
                if '.xml' in hit_info['query_id']:
                    for i in hit_info['query_id'].split('/'):
                        if '.xml' in i:
                            qid = i.split('.')[0] + '.' + i.split('.')[1]
                else:
                    qid = hit_info['query_id']

                _Intermediate["#Query_id"] = qid
                _Intermediate['Match_length'] = hit_info['hit_length']
                _Intermediate['e-value'] = hit_info['e-value']
                _Intermediate['Query_Coverage'] = hit_info['qcov']
                _Intermediate['Hit_Coverage'] = hit_info['scov']
            
                Missing_infor = hmmtop_df.loc[hmmtop_df['Hit_xid'] == hit_xid]
                _Intermediate['Query_Pfam'] = str(hit_info['query_pfams'])
                _Intermediate['Subject_Pfam'] = str(hit_info['subject_pfams'])
                written.append(tcid + '-' + hit_xid)
                _Intermediate["Hit_n_TMS"] = Missing_infor["Hit_n_TMS"].iloc[0] if not Missing_infor.empty else "NA"
                _Intermediate['Query_n_TMS'] = hit_info['Query_n_TMS']
                _Intermediate['Hit_n_TMS'] = hit_info['Hit_n_TMS']
                _Intermediate['TM_Overlap_Score'] = hit_info['TM_Overlap_Score']
            else:
            
                Missing_infor = hmmtop_df.loc[hmmtop_df['Hit_xid'] == hit_xid]

                _Intermediate["Match_length"] = Missing_infor["Match_length"].iloc[0] if not Missing_infor.empty else "NA"
                _Intermediate["Hit_n_TMS"] = Missing_infor["Hit_n_TMS"].iloc[0] if not Missing_infor.empty else "NA"
            
            _Intermediate['isFusion'] = 'NA' 
            _Intermediate["Missing_components"] = str(Output_dict['Missing_proteins'])
            with open(filename, mode="a", encoding='utf-8') as f:
                _Intermediate.to_csv(f, sep='\t', header=False, index=False)
                assign_family(_Intermediate, Output_dict['color'])
        
        

def write_singlecomp(output_dict,Output_df_row):
    color = output_dict[0]
    dictionary = output_dict[1]
    Intermediate = Output_df_row.copy()
    Intermediate['isFusion'] = dictionary[len(dictionary) - 1]
    missing_proteins = "NA"
    Intermediate['Missing_components']= missing_proteins
    Intermediate['Initial_decision'] = 'NA'
    filename=f"{color}.tsv"
    filename = GENOME + '/analysis/' + filename

    filemode='a' if os.path.exists(filename) else 'w' 
    #print(Intermediate)
    with open(filename, mode=filemode, encoding='utf-8') as f:
        Intermediate.to_csv(f, sep='\t', header=filemode=='w', index=False)



'''
def remove_duplicates(input_file, output_file):
    seen = set()  # Track seen rows
    with open(input_file, 'r', newline='') as input_file, \
         open(output_file, 'w', newline='') as output_file:
        reader = csv.reader(input_file, delimiter='\t')
        writer = csv.writer(output_file, delimiter='\t')

        for index, row in enumerate(reader):
            if index == 0:  # Write the first row
                writer.writerow(row)
                continue

            row_tuple = tuple(row)  # Convert row to a hashable tuple

            if row_tuple not in seen:
                seen.add(row_tuple)
                writer.writerow(row)
'''
def custom_sort(row):
    parts = row.split('.')
    num1 = int(parts[0])
    letter = parts[1]
    num2 = int(parts[2])
    num3 = int(parts[3])
    num4 = int(parts[4])
    return num1, letter, num2, num3, num4

def decide_pfam(hit_p, query_p):
    if len(hit_p) == 0 or len(query_p) == 0:
        return False

    for h in hit_p:
        if h in query_p:
            return True

    return False

def perform_pfam():
    yellows = family_assignments['Yellow']
    for family in yellows:
        if family not in family_assignments['Green']:
            subfam = family_assignments['Yellow'][family]
            # looks for missing components to see if they are actually found in the genome
            for sys in subfam:
                system = family_assignments['Yellow'][family][sys]
                num_comps = len(system)
                count = 0
                for row in system:
                    prot_id = row['Hit_tcid'] + '-' + row['Hit_xid']
                    if row['#Query_id'] == 'NA':
                        missing_info[prot_id] = {}
    
    print('Finished yellow pfams')
    reds = family_assignments['Red']
    for family in reds:
        if family not in family_assignments['Green']:
            subfam = family_assignments['Red'][family]
            # looks for missing components to see if they are actually found in the genome
            for sys in subfam:
                system = family_assignments['Red'][family][sys]
                num_comps = len(system)
                count = 0
                for row in system:
                    prot_id = row['Hit_tcid'] + '-' + row['Hit_xid']
                    if row['#Query_id'] == 'NA':
                        missing_info[prot_id] = {}




def find_missing_families():
    yellows = family_assignments['Yellow']
    for family in family_assignments['Yellow']:
        if family not in family_assignments['Green']:
            systems_found = 0
            subfam = family_assignments['Yellow'][family]
            # looks for missing components to see if they are actually found in the genome
            for sys in subfam:
                system = family_assignments['Yellow'][family][sys]
                num_comps = len(system)
                count = 0
                green_queries = []
                green_proteins = []
                for row in system:
                    prot_id = row['Hit_tcid'] + '-' + row['Hit_xid']
                    if prot_id in tcdb_pfams:
                        hit_pfams = tcdb_pfams[prot_id]
                    else:
                        hit_pfams = []
                    if row['#Query_id'] == 'NA':
                        if prot_id not in missing_info:
                            continue
                        hit = missing_info[prot_id]
                        qid = ''

                        if len(hit) == 0:
                            continue
                        for i in hit['query_id'].split('/'):
                            if '.xml' in i:
                                qid = i.split('.')[0] + '.' + i.split('.')[1]
                      
                        # if the protein already has a better match in green continue
                        if row['Hit_xid'] in multi_system_assignments and multi_system_assignments[row['Hit_xid']][1] == 'Green':
                            continue
                        # TODO: fix the duplicate query id error create dict that assigns queries to green and compare w that
                        if qid in queries and queries[qid] == 'Green':
                            continue
                        query_acc = qid
                        if query_acc in genome_pfams:
                            query_pfams = genome_pfams[query_acc]
                        else:
                            query_pfams = []
                        if query_acc not in missing_tms_dict:
                            continue

                        qtms = missing_tms_dict[query_acc]
                        ttms = missing_tms_dict[prot_id]
                        overlap = missing_overlap_dict[query_acc]['alignedTMS']
                        if hit['hit_length'] > MIN_HIT_THRESH and hit['e-value'] <= E_VAL_GREEN and (hit['qcov'] >= Q_COV_THRESH and hit['scov'] >= S_COV_THRESH) and ((qtms > 1 or ttms > 1) and float(overlap / max(qtms, ttms) > 0.5)):
                            multi_system_assignments[row['Hit_xid']][1] = 'Green'
                            count += 1
                            hit['query_pfams'] = query_pfams
                            hit['subject_pfams'] = hit_pfams
                            hit['Query_n_TMS'] = qtms
                            hit['Hit_n_TMS'] = ttms
                            hit['TM_Overlap_Score'] = overlap
                            completed_info[prot_id] = hit
                            green_queries.append(query_acc)
                            green_proteins.append(row['Hit_xid'])
                        elif hit['hit_length'] < MIN_HIT_THRESH and hit['e-value'] <= E_VAL_YELLOW and hit['qcov'] >= LOW_COV and hit['scov'] >= LOW_COV: 
                            multi_system_assignments[row['Hit_xid']][1] = 'Green'
                            count += 1
                            hit['query_pfams'] = query_pfams
                            hit['subject_pfams'] = hit_pfams
                            hit['Query_n_TMS'] = qtms
                            hit['Hit_n_TMS'] = ttms
                            hit['TM_Overlap_Score'] = overlap
                            completed_info[prot_id] = hit
                            green_queries.append(query_acc)
                            green_proteins.append(row['Hit_xid'])
                        elif hit['e-value'] <= E_VAL_YELLOW and decide_pfam(hit_pfams, query_pfams):
                            multi_system_assignments[row['Hit_xid']][1] = 'Green'
                            count += 1
                            hit['query_pfams'] = query_pfams
                            hit['subject_pfams'] = hit_pfams
                            hit['Query_n_TMS'] = qtms
                            hit['Hit_n_TMS'] = ttms
                            hit['TM_Overlap_Score'] = overlap
                            completed_info[prot_id] = hit
                            green_queries.append(query_acc)
                            green_proteins.append(row['Hit_xid'])
                        elif hit['qcov'] >= Q_COV_THRESH and hit['scov'] >= S_COV_THRESH:
                            multi_system_assignments[row['Hit_xid']][1] = 'Green'
                            count += 1
                            hit['query_pfams'] = query_pfams
                            hit['subject_pfams'] = hit_pfams
                            hit['Query_n_TMS'] = qtms
                            hit['Hit_n_TMS'] = ttms
                            hit['TM_Overlap_Score'] = overlap
                            completed_info[prot_id] = hit
                            green_queries.append(query_acc)
                            green_proteins.append(row['Hit_xid'])

                        elif  (qtms > 1 or ttms > 1) and float(overlap / max(qtms, ttms) > 0.5):
                            multi_system_assignments[row['Hit_xid']][1] = 'Green'
                            count += 1
                            hit['query_pfams'] = query_pfams
                            hit['subject_pfams'] = hit_pfams
                            hit['Query_n_TMS'] = qtms
                            hit['Hit_n_TMS'] = ttms
                            hit['TM_Overlap_Score'] = overlap
                            completed_info[prot_id] = hit
                            green_queries.append(query_acc)
                            green_proteins.append(row['Hit_xid'])
                    else:
                        count += 1
                if count == num_comps:
                    systems_found += 1
                    family_assignments['Green'][family] = {}
                    family_assignments['Green'][family][sys] = system
                    family_assignments['Yellow'][family][sys] = None
                    tcid_assignments[sys] = 'Green'
                    move_to_green.append(sys)

                    for i in range(0, len(green_queries)):
                        protein = green_queries[i]
                        queries[protein] = 'Green'
                        matches[protein] = sys
                        tcid_queries[sys].append(protein)
                        curr_hit[sys + '-' + green_proteins[i]]
                else:
                    if family not in family_assignments['Red']:
                        family_assignments['Red'][family] = {}
                    
                    family_assignments['Red'][family][sys] = system
                    family_assignments['Yellow'][family][sys] = None
                    tcid_assignments[sys] = 'Red'
                    move_to_red.append(sys)
               
            if systems_found == 0:
                fam_not_found.append(family)

        # if it is already in greens 
        else:
            subfam = family_assignments['Yellow'][family]
            # looks for missing components to see if they are actually found in the genome
            for sys in subfam:
                system = family_assignments['Yellow'][family][sys]
                if system == None:
                    continue
                num_comps = len(system)
                count = 0
                green_queries = []
                green_proteins = [find_mistake]
                for row in system:
                    prot_id = row['Hit_tcid'] + '-' + row['Hit_xid']
                    if prot_id in tcdb_pfams:
                        hit_pfams = tcdb_pfams[prot_id]
                    else:
                        hit_pfams = []
                    if row['#Query_id'] == 'NA':
                        if prot_id not in missing_info:
                            continue
                        hit = missing_info[prot_id]
                        qid = ''
                        if len(hit) == 0:
                            continue
                        for i in hit['query_id'].split('/'):
                            if '.xml' in i:
                                qid = i.split('.')[0] + '.' + i.split('.')[1]
                        # if the protein already has a better match in green continue
                        if row['Hit_xid'] in multi_system_assignments and multi_system_assignments[row['Hit_xid']][1] == 'Green':
                            continue
                        if qid in queries and queries[qid] == 'Green':
                            continue
                        query_acc = qid
                        if query_acc in genome_pfams:
                            query_pfams = genome_pfams[query_acc]
                        else:
                            query_pfams = []
                        if query_acc not in missing_tms_dict:
                            continue
                        qtms = missing_tms_dict[query_acc]
                        ttms = missing_tms_dict[prot_id]
                        overlap = missing_overlap_dict[query_acc]['alignedTMS']
                        if hit['hit_length'] > MIN_HIT_THRESH and hit['e-value'] <= E_VAL_GREEN and (hit['qcov'] >= Q_COV_THRESH and hit['scov'] >= S_COV_THRESH) and ((qtms > 1 or ttms > 1) and float(overlap / max(qtms, ttms) > 0.5)):
                            hit['query_pfams'] = query_pfams
                            hit['subject_pfams'] = hit_pfams
                            hit['Query_n_TMS'] = qtms
                            hit['Hit_n_TMS'] = ttms
                            hit['TM_Overlap_Score'] = overlap
                            completed_info[prot_id] = hit    
                            count += 1
                            green_queries.append(query_acc)
                            green_proteins.append(row['Hit_xid'])
                        elif hit['hit_length'] < MIN_HIT_THRESH and hit['e-value'] <= E_VAL_YELLOW and hit['qcov'] >= LOW_COV and hit['scov'] >= LOW_COV:
                            count += 1
                            hit['query_pfams'] = query_pfams
                            hit['subject_pfams'] = hit_pfams
                            hit['Query_n_TMS'] = qtms
                            hit['Hit_n_TMS'] = ttms
                            hit['TM_Overlap_Score'] = overlap
                            completed_info[prot_id] = hit
                            green_queries.append(query_acc)
                            green_proteins.append(row['Hit_xid'])
                        elif hit['e-value'] <= E_VAL_YELLOW and decide_pfam(hit_pfams, query_pfams):
                            count += 1
                            hit['query_pfams'] = query_pfams
                            hit['subject_pfams'] = hit_pfams
                            hit['Query_n_TMS'] = qtms
                            hit['Hit_n_TMS'] = ttms
                            hit['TM_Overlap_Score'] = overlap
                            completed_info[prot_id] = hit
                            green_queries.append(query_acc)
                            green_proteins.append(row['Hit_xid'])
                        elif hit['qcov'] >= Q_COV_THRESH and hit['scov'] >= S_COV_THRESH:
                            multi_system_assignments[row['Hit_xid']][1] = 'Green'
                            count += 1
                            hit['query_pfams'] = query_pfams
                            hit['subject_pfams'] = hit_pfams
                            hit['Query_n_TMS'] = qtms
                            hit['Hit_n_TMS'] = ttms
                            hit['TM_Overlap_Score'] = overlap
                            completed_info[prot_id] = hit
                            green_queries.append(query_acc)
                            green_proteins.append(row['Hit_xid'])
                        elif  (qtms > 1 or ttms > 1) and float(overlap / max(qtms, ttms) > 0.5):
                            count += 1
                            hit['query_pfams'] = query_pfams
                            hit['subject_pfams'] = hit_pfams
                            hit['Query_n_TMS'] = qtms
                            hit['Hit_n_TMS'] = ttms
                            hit['TM_Overlap_Score'] = overlap
                            completed_info[prot_id] = hit
                            green_queries.append(query_acc)
                            green_proteins.append(row['Hit_xid'])
                    else:
                        count += 1
                if count == num_comps:
                    family_assignments['Green'][family] = {}
                    family_assignments['Green'][family][sys] = system
                    family_assignments['Yellow'][family][sys] = None
                    tcid_assignments[sys] = 'Green'
                    move_to_green.append(sys)

                    for i in range(0, len(green_queries)):
                        protein = green_queries[i]
                        queries[protein] = 'Green'
                        matches[protein] = sys
                        tcid_queries[sys].append(protein)
                        curr_hit[sys + '-' + green_proteins[i]] = protein
                else:
                    if family not in family_assignments['Red']:
                        family_assignments['Red'][family] = {}
                    
                    family_assignments['Red'][family][sys] = system
                    family_assignments['Yellow'][family][sys] = None
                    tcid_assignments[sys] = 'Red'
                    move_to_red.append(sys)

    reds = family_assignments['Red']
    for family in reds:
        if family not in family_assignments['Green']:
            subfam = family_assignments['Red'][family]
            systems_found = 0
            # looks for missing components to see if they are actually found in the genome
            for sys in subfam:
                system = family_assignments['Red'][family][sys]
                num_comps = len(system)
                count = 0
                green_queries = []
                green_proteins = []
                for row in system:
                    prot_id = row['Hit_tcid'] + '-' + row['Hit_xid']
                    if prot_id in tcdb_pfams:
                        hit_pfams = tcdb_pfams[prot_id]
                    else:
                        hit_pfams = []
                    if row['#Query_id'] == 'NA':
                        if prot_id not in missing_info:
                            continue
                        hit = missing_info[prot_id]
                        qid = ''
                        if len(hit) == 0:
                            continue
                        for i in hit['query_id'].split('/'):
                            if '.xml' in i:
                                qid = i.split('.')[0] + '.' + i.split('.')[1]
                        # if the protein already has a better match in green continue
                        if row['Hit_xid'] in multi_system_assignments and multi_system_assignments[row['Hit_xid']][1] == 'Green':
                            continue
                        if qid in queries and queries[qid] == 'Green' and qid not in matches:
                            continue
                        query_acc = qid
                        if query_acc in genome_pfams:
                            query_pfams = genome_pfams[query_acc]
                        else:
                            query_pfams = []
                        if query_acc in missing_tms_dict:
                            qtms = missing_tms_dict[query_acc]
                        else:
                            qtms = 0
                        if prot_id in missing_tms_dict:
                            ttms = missing_tms_dict[prot_id]
                        else:
                            ttms = 0
                        if query_acc in missing_overlap_dict:
                            overlap = missing_overlap_dict[query_acc]['alignedTMS']
                        else:
                            overlap = 0
                        if hit['hit_length'] > MIN_HIT_THRESH and hit['e-value'] <= E_VAL_GREEN and (hit['qcov'] >= Q_COV_THRESH and hit['scov'] >= S_COV_THRESH):
                            hit['query_pfams'] = query_pfams
                            hit['subject_pfams'] = hit_pfams
                            hit['Query_n_TMS'] = qtms
                            hit['Hit_n_TMS'] = ttms
                            hit['TM_Overlap_Score'] = overlap
                            completed_info[prot_id] = hit
                            count += 1
                            green_queries.append(query_acc)
                            green_proteins.append(row['Hit_xid'])
                        elif hit['hit_length'] < MIN_HIT_THRESH and hit['e-value'] <= E_VAL_YELLOW and hit['qcov'] >= LOW_COV and hit['scov'] >= LOW_COV:
                            count += 1
                            hit['query_pfams'] = query_pfams
                            hit['subject_pfams'] = hit_pfams
                            hit['Query_n_TMS'] = qtms
                            hit['Hit_n_TMS'] = ttms
                            hit['TM_Overlap_Score'] = overlap
                            completed_info[prot_id] = hit
                            green_queries.append(query_acc)
                            green_proteins.append(row['Hit_xid'])
                        elif hit['qcov'] >= Q_COV_THRESH and hit['scov'] >= S_COV_THRESH:
                            count += 1
                            hit['query_pfams'] = query_pfams
                            hit['subject_pfams'] = hit_pfams
                            hit['Query_n_TMS'] = qtms
                            hit['Hit_n_TMS'] = ttms
                            hit['TM_Overlap_Score'] = overlap
                            completed_info[prot_id] = hit
                            green_queries.append(query_acc)
                            green_proteins.append(row['Hit_xid'])
                        elif hit['e-value'] <= E_VAL_YELLOW and decide_pfam(hit_pfams, query_pfams):
                            count += 1
                            hit['query_pfams'] = query_pfams
                            hit['subject_pfams'] = hit_pfams
                            hit['Query_n_TMS'] = qtms
                            hit['Hit_n_TMS'] = ttms
                            hit['TM_Overlap_Score'] = overlap
                            completed_info[prot_id] = hit
                            green_queries.append(query_acc)
                            green_proteins.append(row['Hit_xid'])
                    else:
                        count += 1
                if count == num_comps:
                    systems_found += 1
                    family_assignments['Green'][family] = {}
                    family_assignments['Green'][family][sys] = system
                    family_assignments['Red'][family][sys] = None
                    tcid_assignments[sys] = 'Green'
                    move_to_green.append(sys)

                    for i in range(0, len(green_queries)):
                        protein = green_queries[i]
                        queries[protein] = 'Green'
                        matches[protein] = sys
                        curr_hit[sys + '-' + green_proteins[i]] = protein
            
            if systems_found == 0:
                fam_not_found.append(family)



def find_essential(essentials):
    found = False
    for family in essentials:
        
        # if the family is not already in the greens, it has to be moved there somehow
        if family not in family_assignments['Green']:
            # checks first for the family in the yellows as they are the best hits
            if family in family_assignments['Yellow']:
                subfam = family_assignments['Yellow'][family]
                # looks for missing components to see if they are actually found in the genome
                for sys in subfam:
                    # extract the rows correlating to the system
                    system = family_assignments['Yellow'][family][sys]
                    num_comps = len(system)
                    count = 0
                    green_queries = []
                    green_proteins = []

                    for row in system:
                        prot_id = row['Hit_tcid'] + '-' + row['Hit_xid']
                        # extract the tcdb pfam domains for the specific protein id
                        if prot_id in tcdb_pfams:
                            hit_pfams = tcdb_pfams[prot_id]
                        else:
                            hit_pfams = []
                        # we are looking for specifically missing components to be reclassified so in this case their values are all NA
                        if row['#Query_id'] == 'NA':
                            hit = missing_info[prot_id]
                            qid = ''
                            # extract query id from the xml format
                            for i in hit['query_id'].split('/'):
                                if '.xml' in i:
                                    qid = i.split('.')[0] + '.' + i.split('.')[1]
                            # if the component is already assigned we must complete this system so remove the compoent from the other system
                            if qid in queries and queries[qid] == 'Green':
                                tcdbSystems[row['Hit_tcid']].remove(prot_id)
                                queries[qid] = None
                            
                            query_acc = qid
                            
                            if query_acc not in missing_tms_dict:
                                continue

                            qtms = missing_tms_dict[query_acc]
                            ttms = missing_tms_dict[prot_id]
                            overlap = missing_overlap_dict[query_acc]['alignedTMS']

                            # extract pfams for the query proteins
                            if query_acc in genome_pfams:
                                query_pfams = genome_pfams[query_acc]
                            else:
                                query_pfams = []
                            
                            # if it is a large hit it must have a good e value and good coverage
                            if hit['hit_length'] > MIN_HIT_THRESH and hit['e-value'] <= E_VAL_GREEN and (hit['qcov'] >= Q_COV_THRESH and hit['scov'] >= S_COV_THRESH):
                                count += 1
                                hit['query_pfams'] = query_pfams
                                hit['subject_pfams'] = hit_pfams
                                hit['Query_n_TMS'] = qtms
                                hit['Hit_n_TMS'] = ttms
                                hit['TM_Overlap_Score'] = overlap
                                completed_info[prot_id] = hit
                                green_queries.append(query_acc)
                                green_proteins.append(row['Hit_xid'])
                            # if its a smaller hit it has to have an ok evalue and decent coverage
                            elif hit['hit_length'] < MIN_HIT_THRESH and hit['e-value'] <= E_VAL_YELLOW and hit['qcov'] >= LOW_COV and hit['scov'] >= LOW_COV:
                                hit['query_pfams'] = query_pfams
                                hit['subject_pfams'] = hit_pfams
                                hit['Query_n_TMS'] = qtms
                                hit['Hit_n_TMS'] = ttms
                                hit['TM_Overlap_Score'] = overlap
                                completed_info[prot_id] = hit
                                count += 1
                                green_queries.append(query_acc)
                                green_proteins.append(row['Hit_xid'])
                            elif hit['qcov'] >= Q_COV_THRESH and hit['scov'] >= S_COV_THRESH:
                                count += 1
                                hit['query_pfams'] = query_pfams
                                hit['subject_pfams'] = hit_pfams
                                hit['Query_n_TMS'] = qtms
                                hit['Hit_n_TMS'] = ttms
                                hit['TM_Overlap_Score'] = overlap
                                completed_info[prot_id] = hit
                                green_queries.append(query_acc)
                                green_proteins.append(row['Hit_xid'])
                            # if it has an ok evalue and theres a matching pfam it can be a match
                            elif hit['e-value'] <= E_VAL_YELLOW and decide_pfam(hit_pfams, query_pfams):
                                count += 1
                                hit['query_pfams'] = query_pfams
                                hit['subject_pfams'] = hit_pfams
                                hit['Query_n_TMS'] = qtms
                                hit['Hit_n_TMS'] = ttms
                                hit['TM_Overlap_Score'] = overlap
                                completed_info[prot_id] = hit
                                green_queries.append(query_acc)
                                green_proteins.append(row['Hit_xid'])
                        else:
                            count += 1
                    # reassign component based on if all the components are found
                    if count == num_comps:
                        family_assignments['Green'][family] = {}
                        family_assignments['Green'][family][sys] = system
                        family_assignments['Yellow'][family][sys] = None
                        tcid_assignments[sys] = 'Green'
                        for i in range(0, len(green_queries)):
                            protein = green_queries[i]
                            matches[protein] = sys
                            curr_hit[sys + '-' + green_proteins[i]] = protein
            # if its only found in the red we have to find the system with the most components found and complete that
            else:

                # for the red we need to move the family with the most found components
                subfam = family_assignments['Red'][family]
                least_missing_components = []
                # updates least missing components for each system
                #print(list(subfam.keys()))
                for system in subfam:
                    least_missing_components.append((len(subfam[system][0]['Missing_components']), system))
                sorted_comps = []
                if len(least_missing_components) > 1:
                    sorted_comps = sorted(least_missing_components, key=lambda x:x[0])
                else:
                    sorted_comps = least_missing_components
                #print(sorted_comps)
                for comp in sorted_comps:
                    tcid = comp[1]
                    system = family_assignments['Red'][family][tcid]
                    count = 0
                    num_comps = len(system)
                    green_queries = []
                    green_proteins = []
                    for row in system:
                        prot_id = row['Hit_tcid'] + '-' + row['Hit_xid']
                        if row['#Query_id'] == 'NA':
                            hit = missing_info[prot_id]
                            if len(hit) == 0:
                                continue
                            qid = ''
                            for i in hit['query_id'].split('/'):
                                if '.xml' in i:
                                    qid = i.split('.')[0] + '.' + i.split('.')[1]
 
                            query_acc = qid
                            
                            query_pfams = []
                            hit_pfams = []
                            qtms = 0
                            ttms = 0
                            overlap = 0

                            if query_acc in genome_pfams:
                                query_pfams = genome_pfams[query_acc]
                            if prot_id in tcdb_pfams:
                                hit_pfams = tcdb_pfams[prot_id]
                            if query_acc in missing_tms_dict:
                                qtms = missing_tms_dict[query_acc]
                            if prot_id in missing_tms_dict:
                                ttms = missing_tms_dict[prot_id]
                            if query_acc in missing_overlap_dict:
                                overlap = missing_overlap_dict[query_acc]['alignedTMS']
                            if hit['hit_length'] > MIN_HIT_THRESH and hit['e-value'] <= E_VAL_GREEN and (hit['qcov'] >= Q_COV_THRESH and hit['scov'] >= S_COV_THRESH):
                                count += 1
                                hit['query_pfams'] = query_pfams
                                hit['subject_pfams'] = hit_pfams
                                hit['Query_n_TMS'] = qtms
                                hit['Hit_n_TMS'] = ttms
                                hit['TM_Overlap_Score'] = overlap
                                completed_info[prot_id] = hit
                                green_queries.append(query_acc)
                                green_proteins.append(row['Hit_xid'])
                            elif hit['hit_length'] < MIN_HIT_THRESH and hit['e-value'] <= E_VAL_YELLOW and hit['qcov'] >= LOW_COV and hit['scov'] >= LOW_COV:
                                hit['query_pfams'] = query_pfams
                                hit['subject_pfams'] = hit_pfams
                                hit['Query_n_TMS'] = qtms
                                hit['Hit_n_TMS'] = ttms
                                hit['TM_Overlap_Score'] = overlap
                                completed_info[prot_id] = hit
                                count += 1
                                green_queries.append(query_acc)
                                green_proteins.append(row['Hit_xid'])
                            elif hit['e-value'] <= E_VAL_YELLOW and decide_pfam(hit_pfams, query_pfams):
                                hit['query_pfams'] = query_pfams
                                hit['subject_pfams'] = hit_pfams
                                hit['Query_n_TMS'] = qtms
                                hit['Hit_n_TMS'] = ttms
                                hit['TM_Overlap_Score'] = overlap
                                completed_info[prot_id] = hit
                                count += 1
                                green_proteins.append(row['Hit_xid'])
                                green_queries.append(query_acc)
                            elif hit['qcov'] >= Q_COV_THRESH and hit['scov'] >= S_COV_THRESH:
                                multi_system_assignments[row['Hit_xid']][1] = 'Green'
                                count += 1
                                hit['query_pfams'] = query_pfams
                                hit['subject_pfams'] = hit_pfams
                                hit['Query_n_TMS'] = qtms
                                hit['Hit_n_TMS'] = ttms
                                hit['TM_Overlap_Score'] = overlap
                                completed_info[prot_id] = hit
                                green_queries.append(query_acc)
                                green_proteins.append(row['Hit_xid'])
                        else:
                            count += 1
                    if num_comps - count <= 3:
                        family_assignments['Green'][family] = {}
                        family_assignments['Green'][family][tcid] = system
                        family_assignments['Red'][family][tcid] = None
                        tcid_assignments[tcid] = 'Green'

                        for i in range(0, len(green_queries)):
                            protein = green_queries[i]
                            matches[protein] = tcid
                            curr_hit[tcid + '-' + green_proteins[i]] = protein
                        found = True

    return found

def main():
    # GENOME = 'MicrobiomeResults/GCF_000013425.1'  # we can make this a command line argument
    
    parser = argparse.ArgumentParser(description='Microbiome analysis script')
    parser.add_argument('-g', '--genome', type=str, help='Genome ID')
    parser.add_argument('-q', '--qcov', type=float, help='Query coverage threshold')
    parser.add_argument('-s', '--scov', type=float, help='Subject coverage threshold')
    parser.add_argument('-r', '--autored', type=float, help='Lowest possible coverage allowed')
    parser.add_argument('-m', '--membrane', type=int, help='Membrane proteins threshold')
    parser.add_argument('-t', '--tmsdiff', type=int, help='Lowest amount of TMSs to be a Membrane Protein')
    parser.add_argument('-o', '--tmsthresh', type=float, help='Minimum overlap percentage to be a good hit')
    parser.add_argument('-f', '--feattable', type=str, help='Path to genomic feature table')
    parser.add_argument('-b', '--gbf', type=str, help='Path to gbff file')

    args = parser.parse_args()
    '''
    if not os.path.exists(args.genome):
        print(f"Genome folder not found: {args.genome}")
        exit()
    '''
    global GENOME
    global Q_COV_THRESH
    global S_COV_THRESH
    global AUTO_RED
    global Membraneprotein_threshold
    global Hit_TMS_Diff
    global df
    global TMS_THRESH_GREEN
    global FEATURE_TABLE
    global GBFF_FILE
    if args.genome is not None:
        GENOME = args.genome
    elif len(GENOME) == 0:
        return 'Genome not provided'
    if args.qcov is not None:
        Q_COV_THRESH = args.qcov
    if args.scov is not None:
        S_COV_THRESH = args.scov
    if args.autored is not None:
        AUTO_RED = args.autored
    if args.membrane is not None:
        Membraneprotein_threshold = args.membrane
    if args.tmsdiff is not None:
        Hit_TMS_Diff = args.tmsdiff
    if args.tmsthresh is not None:
        TMS_THRESH_GREEN = args.tmsthresh
    if args.feattable is not None:
        FEATURE_TABLE = args.feattable
    if args.gbf is not None:
        GBFF_FILE = args.gbf
    print('starting')
    # construct df with all data from results.tsv
    df = pd.read_table(GENOME + '/results.tsv')
    setGenome(GENOME)
    # print columns of df for dev use in constructing filtered_df later
    print(GENOME)
    if not os.path.exists(GENOME + '/analysis/'):
        os.mkdir(GENOME + '/analysis')
    
    # create empty dfs for tagging with green, yellow, and red labels
    green_df = df.copy()
    green_df = green_df.iloc[0:0]
    yellow_df = green_df.copy()
    red_df = green_df.copy()

    # construct filtered_df with only relevant data points for each protein matche
    filtered_df = df[['#Query_id', '%_identity', 'e-value', 'Q_start', 'Q_end', 'S_start', 'S_end', 
    'Query_Coverage', 'Hit_Coverage', 'Query_Pfam']]
    print('adjusting overlap score')
    adjustOverlapScore(df)

    Output_df= df[['Hit_tcid','Hit_xid','#Query_id','Match_length','e-value','%_identity','Query_Length','Hit_Length','Q_start',
    'Q_end','S_start','S_end','Query_Coverage','Hit_Coverage','Query_n_TMS','Hit_n_TMS','TM_Overlap_Score','Family_Abrv'
    ,'Predicted_Substrate','Query_Pfam','Subject_Pfam']]

    df.to_csv('adj.csv', index=False)

    # To distinguish between single and multiple component systems. 
    # Example content:
    #  1.A.1.1.1 =>  ["1.A.1.1.1-P0A334"],
    #  3.A.1.1.1 =>  ["3.A.1.1.1-P02916", ""3.A.1.1.1-XXXXX", "3.A.1.1.1-YYYYYY", "3.A.1.1.1-ZZZZZZZZZ"],
    #  ....
    
    pwd = GENOME + '/analysis/tcdb.faa'
    if not os.path.exists(pwd):
        cmd = f"cd {GENOME}; extractTCDB.pl -i tcdb -o analysis/ -f fasta"
        execute_command(cmd)
    input = open(GENOME + '/analysis/tcdb.faa')

    parseTCDBcontent(input)
    assign_tcdb_pfams()

    ##with open('hmmtop.out') as f:
    hmmtop_path=GENOME+ "/analysis/hmmtop.out"
    ##Test db file
    '''
    with open(hmmtop_path_db, "rb") as file:
            
            db_data = pic.load(file)
            print(type(db_data))
            csv_file_path = "hmmtop_test_db.csv"  
            pd.DataFrame.from_dict(db_data).to_csv(csv_file_path, index=False)
    '''
    if not os.path.exists(hmmtop_path):
        TCDB_seqs = GENOME + "/analysis/mcs_tcids"
        if os.path.exists(TCDB_seqs):
            os.system(f"rm -r {TCDB_seqs}")
        os.mkdir(TCDB_seqs)
        data = pd.read_csv(GENOME + '/results.tsv', sep='\t')
        hit_tcid_array = data['Hit_tcid'].unique()
        try:
            with open(GENOME + '/analysis/mcs_tcids.txt', 'w') as f:
                for tcid in hit_tcid_array:
                    f.write(tcid + '\n')
                f.close()
        except Exception as e:
            print("An error occurred:", str(e))
        command_1=[f"extractTCDB.pl -i {GENOME}/analysis/mcs_tcids.txt -o {GENOME}/analysis/mcs_tcids -f fasta"]
        joined_commands = ';'.join(command_1)
        execute_command(joined_commands)
        command2=f"cd {GENOME};cat analysis/mcs_tcids/*faa > analysis/all.faa &&hmmtop -if=analysis/all.faa -of=analysis/hmmtop.out -sf=FAS -pi=spred -is=pseudo"
        execute_command(command2)
        
        
    with open(hmmtop_path) as f:
        lines = f.readlines()
        ##Create a dataframe that include above columns 
        for line in lines:
            fields = re.split(r'[ -]+', line.strip())
            ##split the system and protein names
            new_col=[fields[2],fields[3],fields[5]]
            new_row = pd.Series([fields[2],fields[3],fields[5],fields[1]], index=hmmtop_df.columns)
            hmmtop_df.loc[len(hmmtop_df)] = new_row
        hmmtop_df['Hit_n_TMS'] = hmmtop_df['Hit_n_TMS'].astype(int)
    
    parseTCDBcontent(input)
    Missing_protein_list=[]
    genDict(geneFusions, GENOME)
    for filename in [GENOME + "/analysis/Green.tsv",GENOME + "/analysis/Red.tsv",GENOME + "/analysis/Yellow.tsv"]:
        if os.path.exists(filename):
            os.remove(filename)
            
    print('Beginning file write')
    file_found = False
    if os.path.exists(GENOME + '/analysis/query_ids.txt'):
        file_found = True

    # first iteration 
    for index, row in df.iterrows():
        single = isSingleComp(row)

        # if row['Hit_tcid'] == '1.C.3.4.2':
        #     print('here')
        if(not single):
            # Output_dict = isMultiComp(row, df, 0.5, file_found)
            # assign_multicomp(Output_dict, Output_df.loc[[index], Output_df.columns])
            # Write_multicomp(Output_dict, Output_df.loc[[index], Output_df.columns])
            assign_multicomp_sys(row, df, 0.5, file_found)         
        else:
            output_dict = categorizeSingleComp(row)
            if row['Hit_tcid'] == '1.A.22.1.1':
                print(tcid_assignments['1.A.22.1.1'])
            write_singlecomp(output_dict, Output_df.loc[[index],Output_df.columns])
   
    print('Finished writing to files')
    input_files = ['Green.tsv', 'Red.tsv', 'Yellow.tsv']
    output_files = ['Green_adj.tsv', 'Red_adj.tsv', 'Yellow_adj.tsv']
    
    print('Running genomic context')
    getGenomicContext(10, 4)
    print('Running query pfam')
    #second iteration for MCS
    for index, row in df.iterrows():
        single = isSingleComp(row)

        if not single:
            Output_dict = isMultiCompSecondIteration(row, df, 0.5, file_found)
            assign_multicomp(Output_dict, Output_df.loc[[index], Output_df.columns])
    perform_pfam()
    print('Running pfam program')
    print('created both pfam dicts')
    missing_g_dict, missing_g_overlap_dict, missing_g_tms_dict = find_mistake(GENOME, new_matches, 8)
    missing_dict, missing_overlap_dict, missing_tms_dict = find(GENOME, list(missing_info.keys()), 8)

    for m in missing_g_dict:

        if len(missing_g_dict[m]) == 0:
            continue
        
        query = missing_g_dict[m]
        query_id = query['query_id']
        query_pfams = []
        if query_id in genome_pfams:
            query_pfams = genome_pfams[query_id]

        hit_pfams = []
        if m in tcdb_pfams:
            hit_pfams = tcdb_pfams[m]

        missing_info[m] = query
        queries[query['query_id']] = 'Green'
        query['query_pfams'] = query_pfams
        query['subject_pfams'] = hit_pfams
        if query_id in missing_g_tms_dict:
            query['Query_n_TMS'] = missing_g_tms_dict[query_id]
            
        else:
            query['Query_n_TMS'] = 0
        if m in missing_g_tms_dict:
            query['Hit_n_TMS'] = missing_g_tms_dict[m]
        else:
            query['Hit_n_TMS'] = 0

        if query_id in missing_g_overlap_dict:
            query['TM_Overlap_Score'] = missing_g_overlap_dict[query_id]['alignedTMS']
        else:
            query['TM_Overlap_Score'] = 0

        completed_info[m] = query
    for m in missing_dict.keys():
        found_key = False
        count = 0
        if len(missing_dict[m]) == 0:
            continue
        while not found_key:
            query = missing_dict[m][count][0]

            if query['query_id'] not in queries or queries[query['query_id']] != 'Green':
                missing_info[m] = query
                queries[query['query_id']] = 'Green'
                found_key = True
            
            count += 1

            if count == len(missing_dict[m]):
                found_key = True

    
    print('Running pfam program')
    run_query_pfam(file_found)
    
    print('finding essential systems')
    find_essential(['3.A.2', '3.A.5'])
    print('looking for fams')

    find_missing_families()


    print('beginning second iteration')
    
    with open(GENOME + '/analysis/interesting_cases.txt', 'w') as f:
        f.write('Families not present in genome\n')
        if len(fam_not_found) == 0:
            f.write('None\n')
        for fam in fam_not_found:
            f.write(fam + '\n')

        f.write('Systems moved to Green\n')
        for sys in move_to_green:
            f.write(sys + '\n')

        f.write('Systems moved to Red\n')
        for sys in move_to_red:
            f.write(sys + '\n')


    # second iteration
    for index, row in df.iterrows():
        if not isSingleComp(row):
            Output_dict = isMultiComp(row, df, 0.5, file_found)
            Write_multicomp(Output_dict, Output_df.loc[[index], Output_df.columns])
        else:
            output_dict = categorizeSingleComp(row)
    
    for filename in [GENOME + "/analysis/Green.tsv",GENOME + "/analysis/Red.tsv",GENOME + "/analysis/Yellow.tsv"]: 
        if os.path.exists(filename):
            df = pd.read_csv(filename, sep='\t')
            df = df.fillna('NA')
            NA_count = df.eq('NA').sum(axis=1).rename('NA_count')
            df['na_sort'] = NA_count
            df = df.sort_values(by=['Hit_tcid','na_sort']).drop(columns=['na_sort'])
            mask = (df['#Query_id'] == 'NA') & (df['e-value'] == 'NA')
            df_temp = df.loc[mask] 
            df_temp.drop_duplicates(subset=['Hit_tcid','Hit_xid'], keep='first', inplace=True) 
            df.loc[mask] = df_temp
            df.dropna(inplace=True)
            df.to_csv(filename, sep='\t', index=False)
            # df.sort_values(by=df.columns[0], key=lambda x: x.map(custom_sort))
            command = f"sort -o {filename} -t '.' -k1,1n -k2,2 -k3,3n -k4,4n -k5,5n {filename}"
            subprocess.run(command, shell=True)

main()
