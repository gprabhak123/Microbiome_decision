from genome_comparison import *
import pandas as pd
import os
import csv

genomes = []

# generates the columns needed for our Master Dataframe
def master_df_column_gen(directory):
    columns = ['#TCID', 'Acc', 'CE', 'Role', 'hit_tms_no']
    per_genome_columns = ['query', 'q_tms', 'evalue', 'pident', 'qcov', 'scov']
    num_genomes = 0
    items = os.listdir(directory)
    
    global genomes
    genomes = [item for item in items if (os.path.isdir(os.path.join(directory, item)) and 'GCF' in item)]
    
    columns.extend(genomes)
    for genome in genomes:
        columns.extend(per_genome_columns)

    return columns

# declare master_dict as global
master_dict = {}

def master_dict_generation(green_dir):
    getSmithWaterman('/Users/gautham/microbiome_project/test_comparisons')
    parse_sw('/Users/gautham/microbiome_project/test_comparisons')
    for filename in os.listdir(green_dir):
        if os.path.isfile(os.path.join(green_dir, filename)):

            file_parts = filename.split('_')[:2]
            genome = '_'.join(file_parts)

            green_df = pd.read_csv(os.path.join(green_dir, filename), delimiter='\t')
            
            for index, row in green_df.iterrows():
                tcid_acc = row['Hit_tcid'] + '-' + row['Hit_xid']

                try:
                    info_df = getInfoAsRow(genome, tcid_acc)
                except:
                    fix_cmd = f"grep {row['Hit_xid']} ~/db/blastdb/tcdb.faa"
                    correct_tcdb = subprocess.run([fix_cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True).stdout.strip().split('>')[1]
                    tcid_acc = correct_tcdb
                    info_df = getInfoAsRow(genome, tcid_acc)

                if tcid_acc in master_dict and master_dict[tcid_acc]:
                    if genome in master_dict[tcid_acc] and master_dict[tcid_acc][genome]:
                        master_dict[tcid_acc][genome][row['#Query_id']] = {}
                    else:
                        master_dict[tcid_acc][genome] = {row['#Query_id']:{}}
                else: 
                    master_dict[tcid_acc] = {genome: {row['#Query_id']:{}}}
                
                info_row = info_df.iloc[0]
               
                qcov_val = round(float(((info_row['hit_end'] - info_row['hit_start'] + 1) / info_row['hit_len']) * 100), 1)
                scov_val = round(float(((info_row['query_end'] - info_row['query_start'] + 1) / info_row['query_len']) * 100), 1)

                #print(f"QID: {row['#Query_id']}")
                #print(f"QUERY {qcov_val} - {info_row['hit_end']} - {info_row['hit_start']} + 1) / {info_row['hit_len']}")
                #print(f"SUBJECT {scov_val}  - {info_row['query_end']} - {info_row['query_start']} + 1) / {info_row['query_len']}")
                #qcov_val = int(info_row['hit_end']) - int(info_row['hit_start'])
                #scov_val = int(info_row['query_end']) - int(info_row['query_start'])

                

                master_dict[tcid_acc][genome][row['#Query_id']] = {
                    'CE' : 'N/A',
                    'Role' : 'N/A',
                    'hit_tms_no' : row['Hit_n_TMS'],
                    'query': info_row['hit_id'], 
                    'q_tms': row['Query_n_TMS'],
                    'evalue' :info_row['eval'],
                    'pident':info_row['pident'],
                    'qcov': qcov_val,
                    'scov': scov_val,
                }

NUM_COL_PRE_GENOME = 5
NUM_COL_PER_GENOME = 6

def tsv_generation():
    df_columns = master_df_column_gen('test_comparisons')
    master_dict_generation('test_comparisons/greens')
    master_df = pd.DataFrame(columns=df_columns)

    master_array = []
    master_array.append(df_columns)

    for tcid_acc_key in master_dict:
        output_row = [None] * len(df_columns)
        output_row[0] = tcid_acc_key.split('-')[0]
        output_row[1] = tcid_acc_key.split('-')[1]

        for i in range(len(genomes)):
            if genomes[i] in master_dict[tcid_acc_key]:
                output_row[NUM_COL_PRE_GENOME+i] = '+'
                first_protein = list(master_dict[tcid_acc_key][genomes[i]])[0]
                first_protein_vals = master_dict[tcid_acc_key][genomes[i]][first_protein]
                start_index = NUM_COL_PRE_GENOME+len(genomes)+NUM_COL_PER_GENOME*i-1
                output_row[start_index+1] = first_protein_vals['query']
                output_row[start_index+2] = first_protein_vals['q_tms']
                output_row[start_index+3] = first_protein_vals['evalue']
                output_row[start_index+4] = first_protein_vals['pident']
                output_row[start_index+5] = first_protein_vals['qcov']
                output_row[start_index+6] = first_protein_vals['scov']
                output_row[2] = first_protein_vals['CE']
                output_row[3] = first_protein_vals['Role']
            else:
                output_row[NUM_COL_PRE_GENOME+i] = '-'
                start_index = NUM_COL_PRE_GENOME+len(genomes)+NUM_COL_PER_GENOME*i
                for j in range(6):
                    output_row[start_index+j] = 'none'

        master_array.append(output_row)

    with open('test_array_output.tsv', 'w', newline='') as tsv_output:
        writer = csv.writer(tsv_output, delimiter='\t')
        writer.writerows(master_array)

tsv_generation()
