import xml.etree.ElementTree as ET
import xml
import pandas as pd
import pickle as pic
import json
import subprocess
import os
from mmseqs_hmmtop import overlapDict
import hmmtop


query_data = {}
target_data = {}
overlap_dict = {}

def read_pickle(hmmtop_file):
    global query_data
    global target_data

    with open(hmmtop_file, "rb") as file:
        # Deserialize the data using pickle.load()
        data = pic.load(file)

    # Extracts all query data and target data from hmmtop file
    query_data = data['queries']

    with open('genome_tms.txt', 'w') as f:
        json.dump(data, f, indent=4)

    for key in data['tcdb']:
        key_elems = key.split('|')
        new_key = key_elems[3] + '-' + key_elems[2].split('.')[0]

        target_data[new_key] = data['tcdb'][key]
    
def extract_all(genome_path, multicomp_proteins):
    xml_files = []
    for protein in multicomp_proteins:
        xml_files.append(genome_path + '/xml/' + protein + '.xml')
    
    all_matches = {}
    for xml_path in xml_files:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        
        qid = ''
        for i in xml_path.split('/'):
            if '.xml' in i:
                qid = i.split('.')[0] + '.' + i.split('.')[1]
        q_len = int(root.find('./BlastOutput_iterations/Iteration/Iteration_query-len').text)
        for item in root.findall('./BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
            hit_info = item.find('Hit_def').text.split('|')
            hit_num = int(item.find('Hit_num').text)
            hit_length = int(item.find('Hit_len').text)
            good_hit = True
            for h_item in item.findall('Hit_hsps/Hsp'):
                e_val = float(h_item.find('Hsp_evalue').text)
                qstart = int(h_item.find('Hsp_query-from').text)
                qend = int(h_item.find('Hsp_query-to').text)
                hstart = int(h_item.find('Hsp_hit-from').text)
                hend = int(h_item.find('Hsp_hit-to').text)
 
                q_cov = round(float(((qend - qstart + 1) / q_len) * 100), 1)
                s_cov = round(float(((hend - hstart + 1) /hit_length) * 100), 1)
 
                query_seq = h_item.find('Hsp_qseq').text
                subject_seq = h_item.find('Hsp_hseq').text
                
                if e_val >= 1e-3 and q_cov < 70 and s_cov < 70:
                    good_hit = False
                    break
            if qid not in all_matches:
                all_matches[qid] = {}
            if good_hit == False:
                break
            

            all_matches[qid][hit_num] = hit_info[3].split(' ')[0] + '-' + hit_info[2]

    return all_matches
            

def find_mistake(genome_path, mistakes, minRes):
    missing_dict = {}
    for p in mistakes:
        protein_dict = []
        query = mistakes[p]
        tcid = p.split('-')[0]
        protein = p.split('-')[1]
        file = genome_path + '/xml/' + query + '.xml'
        
        tree = ET.parse(file)
        root = tree.getroot()

        mmseqs = {}
        hmmtop_dict = {}
        to_run_query = []
        to_run_tcdb = []
        tms_dict = {}
        q_len = int(root.find('./BlastOutput_iterations/Iteration/Iteration_query-len').text)
        for item in root.findall('./BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
            hit_info = item.find('Hit_def').text.split('|')
            hit_num = int(item.find('Hit_num').text)
            hit_length = int(item.find('Hit_len').text)
            if hit_info[2].split('.')[0] == protein and hit_info[3].split(' ')[0] == tcid:
                for h_item in item.findall('Hit_hsps/Hsp'):
                    e_val = float(h_item.find('Hsp_evalue').text)
                    qstart = int(h_item.find('Hsp_query-from').text)
                    qend = int(h_item.find('Hsp_query-to').text)
                    hstart = int(h_item.find('Hsp_hit-from').text)
                    hend = int(h_item.find('Hsp_hit-to').text)

                    q_cov = round(float(((qend - qstart + 1) / q_len) * 100), 1)
                    s_cov = round(float(((hend - hstart + 1) /hit_length) * 100), 1)


                    query_seq = h_item.find('Hsp_qseq').text
                    subject_seq = h_item.find('Hsp_hseq').text

                    qid = query

                    if qid not in mmseqs:
                        mmseqs[qid] = {'qaln': query_seq, 'taln': subject_seq, 'target': tcid + '-' + protein, 'qstart': qstart, 'qend': qend, 'tstart': hstart, 'tend': hend}
                    qtms = {}

                    # if the query tmss already are in the pickle as it mauy be a protein matched with another system, just extract info from there
                    if qid in query_data:
                        qtms['tms'] = list(query_data[qid].values())


                        tms_dict[qid] = len(qtms['tms'])
                        if qid not in hmmtop_dict:
                            hmmtop_dict[qid] = qtms
                    # if not we need to get sequences to run hmmtop later
                    else:
                        to_run_query.append(qid)

                    ttms = {}
                    if (tcid + '-' + protein) in target_data:
                        ttms['tms'] = list(target_data[tcid + '-' + protein].values())
                        tms_dict[tcid + '-' + protein] = len(ttms['tms'])
                        if tcid + '-' + protein not in hmmtop_dict:
                            hmmtop_dict[tcid + '-' + protein] = ttms
                    else:
                        to_run_tcdb.append(tcid + '-' + protein)
                    
                    missing_dict[tcid + '-' + protein] = {'tcdb_protein': protein, 'e-value': e_val, 'qcov': q_cov, 'scov': s_cov, 'hit_length': hit_length, 'query_length': q_len, 'query_id': query}
    
    with open(genome_path + '/analysis/g_accessions_tms_query.txt', 'w') as f:
        for a in to_run_query:
            f.write(a + '\n')
 
    with open(genome_path + '/analysis/g_accessions_tms_tcdb.txt', 'w') as f:
        for a in to_run_tcdb:
            f.write(a + '\n')
   
    cmd1 = f'getseqs nr {genome_path}/analysis/accessions_tms_query.txt > {genome_path}/analysis/all_missing_query.faa'
    print(cmd1)
    cmd2 = f'getseqs tcdb {genome_path}/analysis/accessions_tms_tcdb.txt > {genome_path}/analysis/all_missing_tcdb.faa'
    subprocess.run(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(cmd2)
    subprocess.run(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
 
    ht = hmmtop.tools()
 
    res_q,symbols_q = ht.scan_file(genome_path + '/analysis/all_missing_query.faa')
    res_t, symbols_t = ht.scan_file(genome_path + '/analysis/all_missing_tcdb.faa')
    
    qtms = {}
    for i in range(0, len(symbols_q)):
        if symbols_q[i] in hmmtop_dict:
            continue
        else:
            qtms['tms'] = list(res_q[i].values())
            hmmtop_dict[symbols_q[i]] = qtms
            tms_dict[symbols_q[i]] = len(qtms['tms'])
 
    
    ttms = {}
    for i in range(0, len(symbols_t)):
        if symbols_t[i] not in hmmtop_dict:
            ttms['tms'] = list(res_t[i].values())
            hmmtop_dict[symbols_t[i]] = ttms
            tms_dict[symbols_t[i]] = len(ttms['tms'])
 
    for protein in to_run_query:
        if protein not in symbols_q:
            qtms['tms'] = []
            hmmtop_dict[protein] = qtms
            tms_dict[protein] = 0
 
    for protein in to_run_tcdb:
        if protein not in symbols_t:
            qtms['tms'] = []
            hmmtop_dict[protein] = qtms
            tms_dict[protein] = 0
 
    overlap_dict = overlapDict(mmseqs, hmmtop_dict, minRes)
    return (missing_dict, overlap_dict, tms_dict)


def find(genome_path, missing_comps, minRes):
    with open(genome_path + '/analysis/missing_comps.txt', 'w') as f:
        for protein in missing_comps:
            f.write(protein.split('-')[1] + '\n')
            #f.write(protein + '\n')
    
    cmd = f'grep -f {genome_path}/analysis/missing_comps.txt {genome_path}/xml/*.xml'
    #xml_file = os.popen(cmd).read()
    xml_file = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout
    found_proteins = {}
    temp = xml_file.split('\n')

    hmmtop_file = genome_path + '/hmmtop.db'
    read_pickle(hmmtop_file)
    missing_dict = {}
    mmseqs = {}
    hmmtop_dict = {}
    to_run_query = []
    to_run_tcdb = []
    tms_dict = {}
    

    for f in temp:
        protein_dict = []

        
        if f == '':
            continue
        xml_path = f.split(' ')[0][:-1]
        
        tcid = f.split(' ')[2].split('|')[3]
        protein = f.split(' ')[2].split('|')[2]
        tree = ET.parse(xml_path)
        root = tree.getroot()
    
        q_len = int(root.find('./BlastOutput_iterations/Iteration/Iteration_query-len').text)
        for item in root.findall('./BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
            hit_info = item.find('Hit_def').text.split('|')
            hit_num = int(item.find('Hit_num').text)
            hit_length = int(item.find('Hit_len').text)
            if hit_info[2].split('.')[0] == protein and hit_info[3].split(' ')[0] == tcid:
                for h_item in item.findall('Hit_hsps/Hsp'):
                    e_val = float(h_item.find('Hsp_evalue').text)
                    qstart = int(h_item.find('Hsp_query-from').text)
                    qend = int(h_item.find('Hsp_query-to').text)
                    hstart = int(h_item.find('Hsp_hit-from').text)
                    hend = int(h_item.find('Hsp_hit-to').text)
            
                    q_cov = round(float(((qend - qstart + 1) / q_len) * 100), 1)
                    s_cov = round(float(((hend - hstart + 1) /hit_length) * 100), 1)

                    
                    query_seq = h_item.find('Hsp_qseq').text
                    subject_seq = h_item.find('Hsp_hseq').text
                    
                    qid = ''
                    for i in xml_path.split('/'):
                        if '.xml' in i:
                            qid = i.split('.')[0] + '.' + i.split('.')[1]
                    
                    if qid not in mmseqs:
                        mmseqs[qid] = {'qaln': query_seq, 'taln': subject_seq, 'target': tcid + '-' + protein, 'qstart': qstart, 'qend': qend, 'tstart': hstart, 'tend': hend}
                    qtms = {}

                    # if the query tmss already are in the pickle as it mauy be a protein matched with another system, just extract info from there
                    if qid in query_data:
                        qtms['tms'] = list(query_data[qid].values())

                        
                        tms_dict[qid] = len(qtms['tms'])
                        if qid not in hmmtop_dict:
                            hmmtop_dict[qid] = qtms
                    # if not we need to get sequences to run hmmtop later
                    else:
                        to_run_query.append(qid)
                        
                    ttms = {}
                    if (tcid + '-' + protein) in target_data:
                        ttms['tms'] = list(target_data[tcid + '-' + protein].values())
                        tms_dict[tcid + '-' + protein] = len(ttms['tms'])
                        if tcid + '-' + protein not in hmmtop_dict:
                            hmmtop_dict[tcid + '-' + protein] = ttms
                    else:
                        to_run_tcdb.append(tcid + '-' + protein)

        


                    protein_dict.append(({'query_id': xml_path,'e-value': e_val, 'qcov': q_cov, 'scov': s_cov, 'hit_length': hit_length, 'query_length': q_len, 'protein_match': tcid + '-' + protein}, hit_num))
        sorted_proteins = sorted(protein_dict, key=lambda x:x[1])
        missing_dict[tcid + '-' + protein] = sorted_proteins
    with open(genome_path + '/analysis/accessions_tms_query.txt', 'w') as f:
        for a in to_run_query:
            f.write(a + '\n')

    with open(genome_path + '/analysis/accessions_tms_tcdb.txt', 'w') as f:
        for a in to_run_tcdb:
            f.write(a + '\n')
   
    cmd1 = f'getseqs nr {genome_path}/analysis/accessions_tms_query.txt > {genome_path}/analysis/all_missing_query.faa'
    print(cmd1)
    cmd2 = f'getseqs tcdb {genome_path}/analysis/accessions_tms_tcdb.txt > {genome_path}/analysis/all_missing_tcdb.faa'
    subprocess.run(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(cmd2)
    subprocess.run(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)



    ht = hmmtop.tools()

    res_q,symbols_q = ht.scan_file(genome_path + '/analysis/all_missing_query.faa')
    res_t, symbols_t = ht.scan_file(genome_path + '/analysis/all_missing_tcdb.faa')
    
    qtms = {}
    for i in range(0, len(symbols_q)):
        if symbols_q[i] in hmmtop_dict:
            continue
        else:
            qtms['tms'] = list(res_q[i].values())
            hmmtop_dict[symbols_q[i]] = qtms
            tms_dict[symbols_q[i]] = len(qtms['tms'])

    
    ttms = {}
    for i in range(0, len(symbols_t)):
        if symbols_t[i] not in hmmtop_dict:
            ttms['tms'] = list(res_t[i].values())
            hmmtop_dict[symbols_t[i]] = ttms
            tms_dict[symbols_t[i]] = len(ttms['tms'])

    for protein in to_run_query:
        if protein not in symbols_q:
            qtms['tms'] = []
            hmmtop_dict[protein] = qtms
            tms_dict[protein] = 0

    for protein in to_run_tcdb:
        if protein not in symbols_t:
            qtms['tms'] = []
            hmmtop_dict[protein] = qtms
            tms_dict[protein] = 0

    overlap_dict = overlapDict(mmseqs, hmmtop_dict, minRes)
    return (missing_dict, overlap_dict, tms_dict)


def main():
    missing_comps = []
    with open('/ResearchData/Microbiome/gblast/GCF_001689125.2/analysis/missing_comps.txt', 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        protein = line.strip()
        missing_comps.append(protein)

    missing_dict, overlap_dict, tms_dict = find('/ResearchData/Microbiome/gblast/GCF_001689125.2', missing_comps, 8)

#main()
