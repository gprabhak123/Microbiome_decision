import pandas as pd
import numpy as np
import os

GENOME = ''
def setGenome(genome):
    GENOME = genome


geneFusions = {}
test_fus =  [{'query': 'A', 'qcov': 69.6, 'sstart': 295, 'send': 713, 'scov': 58.2}, {'query': 'B', 'qcov': 65.4, 'sstart': 405, 'send': 709, 'scov': 42.3}, {'query': 'C', 'qcov': 41.4, 'sstart': 481, 'send': 717, 'scov': 32.9}]

# compare overlap between individual fusion candidates to ensure there is not over a 20% overlap 
# (if an individual candidate overlaps with all others with greater than thresh, reject)
MAX_FUSION_COV_OVERLAP = 20
MIN_QUERY_COV_THRESHOLD = 80
MIN_SUBJECT_COV_THRESHOLD = 5
MAX_SUBJECT_COV_THRESHOLD = 90
MIN_FUSION_COVERAGE = 50

# generates a dictionary from the dataframe
def genDict(dict, genome):
    df = pd.read_table(genome + '/results.tsv')
    for index, row in df.iterrows():
        new_entry = {"query": row["#Query_id"], "qcov": row["Query_Coverage"], "sstart": row["S_start"], "send": row["S_end"], "scov": row["Hit_Coverage"], "hit_length": row["Hit_Length"], "tms": row['Query_n_TMS']}
        tcid = row["Hit_tcid"] + "-" + row["Hit_xid"]
        if tcid in dict:
            dict.get(tcid).append(new_entry)
            #print("Appending : " + tcid)
        else:
            dict[tcid] = [new_entry]
            #print("New entry")

def check_overlap(dict1, dict2):
    return (dict1['sstart'] <= dict2['send'] and dict1['send'] >= dict2['sstart']) or (dict1['send'] >= dict2['sstart'] and dict1['sstart'] <= dict2['send'])

# takes in sorted array (subset of dictionary) and 
# outputs list of dictionary (each dictionary is a fusion candidate)
def isFusion(sortedArr):

    
    # output list
    fus_list = []
    overlap_length = 0
    
    # iteratate through each fusion candidate 
    for i in range(len(sortedArr)):
        if sortedArr[i]['scov'] > MAX_SUBJECT_COV_THRESHOLD:
            continue
        elif sortedArr[i]['scov'] < MIN_SUBJECT_COV_THRESHOLD:
            continue
        elif sortedArr[i]['qcov'] < MIN_QUERY_COV_THRESHOLD:
            continue
        else:
            fus_list.append(sortedArr[i])
    tot_length = 0
    overlap_length = 0

    fusion_groups = []
    sorted_fus_list = sorted(fus_list, key=lambda d: d['sstart'])     
    for fus in sorted_fus_list:
        overlap_fus = None         
        for fus_group in fusion_groups:             
            if any(check_overlap(fus, d) for d in fus_group):                 
                overlap_fus = fus                 
                if overlap_fus:                     
                    fus_group.append(overlap_fus)                 
                break         
        if overlap_fus == None:             
            fusion_groups.append([fus])


    #checks for extra fusions if proteins are lined up end to end
    for i in range(len(fus_list)):
        tot_length += (fus_list[i]['send'] - fus_list[i]['sstart'])
        if i < len(fus_list)-1:
            overlap_length += min(fus_list[i+1]['send'], fus_list[i]['send']) - max(fus_list[i+1]['sstart'], fus_list[i]['sstart'])
    tot_length -= overlap_length
    
    if len(fus_list) == 0:
        return []
    elif len(fus_list) < 2:
        return []
    elif (tot_length/fus_list[0]['hit_length'])*100 < MIN_FUSION_COVERAGE:
        return []
    elif len(fusion_groups) > 1:
        return fus_list
    else:
        return []

# PRE-CONDITION: There is good coverage and good evalue already for the fusion candidates
def fusion_TMS_count(fus_list):

    net_TMS = 0
    
    # construct "fusion groups"
    fusion_groups = []
    sorted_fus_list = sorted(fus_list, key=lambda d: d['sstart'])

    for fus in sorted_fus_list:
        overlap_fus = None
        for fus_group in fusion_groups:
            #print(fus_group)
            if any(check_overlap(fus, d) for d in fus_group):
                overlap_fus = fus
                if overlap_fus:
                    fus_group.append(overlap_fus)
                break
        if overlap_fus == None:
            fusion_groups.append([fus])
     

            
    
    # print(sorted_fus_list)
    # iterate over "fusion groups" and return highest tms count
    for fus_group in fusion_groups:
        sorted_fus_group = sorted(fus_group, key=lambda d: d['tms'])
        #print(sorted_fus_group)
        #print(sorted_fus_group)
        #print(type(sorted_fus_group[count])
        net_TMS += sorted_fus_group[0]["tms"]
    
    return net_TMS



def main():
    genDict(geneFusions, GENOME)


    #sorts the genomeFusions dictionary by the genome start index
    num_out = 0
    num_in = 0
    for id in geneFusions:
        if(len(geneFusions[id]) == 1):
            continue
        sortedGeneArr = sorted(geneFusions[id], key=lambda x: x['sstart'])
        # print(sortedGeneArr)
        # x = input()
        num_in+=1
        if len(isFusion(sortedGeneArr)) != 0:
            # print("SUCCESS")
            num_out+=1
    #y = input()
# print("-------")
# print(str(num_out) + " fusions found out of " + str(num_in))
