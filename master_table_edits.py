from genome_comparison import *
import pandas as pd
import os
import csv
import copy

genomes = []
pwd = '/Users/gautham/microbiome_project/test_comparisons'


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
    getSmithWaterman(pwd)
    parse_sw(pwd)
    #getSmithWaterman('/Users/gautham/microbiome_project/test_shit')
    #parse_sw('/Users/gautham/microbiome_project/test_shit')

    
    for filename in os.listdir(green_dir):
        if os.path.isfile(os.path.join(green_dir, filename)):

            file_parts = filename.split('_')[:2]
            genome = '_'.join(file_parts)

            green_df = pd.read_csv(os.path.join(green_dir, filename), delimiter='\t')
            
            for index, row in green_df.iterrows():
                tcid_acc = row['Hit_tcid'] + '-' + row['Hit_xid']

                info_df = None
                try:
                    info_df = getInfoAsRow(genome, tcid_acc)
                except:
                    fix_cmd = f"grep {row['Hit_xid']} ~/db/blastdb/tcdb.faa"
                    correct_tcdb = subprocess.run([fix_cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True).stdout.strip().split('>')[1]
                    tcid_acc = correct_tcdb

                if tcid_acc in master_dict and master_dict[tcid_acc]:
                    if genome in master_dict[tcid_acc] and master_dict[tcid_acc][genome]:
                        master_dict[tcid_acc][genome][row['#Query_id']] = {}
                    else:
                        master_dict[tcid_acc][genome] = {row['#Query_id']:{}}
                else: 
                    master_dict[tcid_acc] = {genome: {row['#Query_id']:{}}}
                if info_df is None:
                    master_dict[tcid_acc][genome][row['#Query_id']] = {
                        'CE' : 'N/A',
                        'Role' : 'N/A',
                        'hit_tms_no' : row['Hit_n_TMS'],
                        'query': 'none',
                        'q_tms': row['Query_n_TMS'],
                        'evalue' : 'none',
                        'pident': 'none',
                        'qcov': 'none',
                        'scov': 'none',
                    } 
                else:

                    def which_row(df, query_to_find):
                        for i, r in df.iterrows():
                            if r['hit_id'] == query_to_find:
                                return i
                        
                        return 0

                    idx = which_row(info_df, row['#Query_id'])
                    info_row = info_df.iloc[idx]
               
                    qcov_val = round(float(((info_row['hit_end'] - info_row['hit_start'] + 1) / info_row['hit_len']) * 100), 1)
                    scov_val = round(float(((info_row['query_end'] - info_row['query_start'] + 1) / info_row['query_len']) * 100), 1)

                    #print(f"QID: {row['#Query_id']}")
                    #print(f"QUERY {qcov_val} - {info_row['hit_end']} - {info_row['hit_start']} + 1) / {info_row['hit_len']}")
                    #print(f"SUBJECT {scov_val}  - {info_row['query_end']} - {info_row['query_start']} + 1) / {info_row['query_len']}")
                    #qcov_val = int(info_row['hit_end']) - int(info_row['hit_start'])
                    #scov_val = int(info_row['query_end']) - int(info_row['query_start'])

                    # if tcid_acc == '1.A.115.1.5-WP_062382148':
                    #     print('here')
                    #     print(row['#Query_id'])


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

def master_tsv_generation():
    df_columns = master_df_column_gen('test_comparisons')
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
                output_row[4] = first_protein_vals['hit_tms_no']
            else:
                output_row[NUM_COL_PRE_GENOME+i] = '-'
                start_index = NUM_COL_PRE_GENOME+len(genomes)+NUM_COL_PER_GENOME*i
                for j in range(6):
                    output_row[start_index+j] = 'none'

        master_array.append(output_row)

    with open('test_array_output.tsv', 'w', newline='') as tsv_output:
        writer = csv.writer(tsv_output, delimiter='\t')
        writer.writerows(master_array)

def additional_row(prev_row, protein_vals):
        additional_row = prev_row
        additional_row[2] = protein_vals['CE']
        additional_row[3] = protein_vals['Role']
        additional_row[4] = protein_vals['hit_tms_no']
        additional_row[5] = protein_vals['query']
        additional_row[6] = protein_vals['q_tms']
        additional_row[7] = protein_vals['evalue']
        additional_row[8] = protein_vals['pident']
        additional_row[9] = protein_vals['qcov']
        additional_row[10] = protein_vals['scov']
        return additional_row

def genome_tsv_generation():
    master_array = []
    genome_columns = ['#TCID','Acc','CE','Role','hit_tms_no','query','q_tms','evalue','pident','qcov','scov']
    master_array.append(genome_columns)

    for genome in genomes:

        genome_master_array = []

        for tcid_acc_key in master_dict:
            if genome in master_dict[tcid_acc_key]:
                output_row = [None] * len(genome_columns)
                output_row[0] = tcid_acc_key.split('-')[0]
                output_row[1] = tcid_acc_key.split('-')[1]
                # if 'WP_062382148' in tcid_acc_key:
                #     # print(getInfoAsRow('GCF_001558775.1', '1.A.115.1.5-WP_062382148'))
                #     # print('-------------')
                #     # print(master_dict[tcid_acc_key][genome])

                #### TODO: Error here
                temp_arr = []
                for protein in master_dict[tcid_acc_key][genome]:
                    temp = copy.copy(output_row)
                    genome_master_array.append(additional_row(temp, master_dict[tcid_acc_key][genome][protein]))
                    #genome_master_array.append(row_to_add)
        
        filename = f"{pwd}/{genome}/test_genome_{genome}.tsv"
        with open(filename, 'w', newline='') as tsv_output:
            writer = csv.writer(tsv_output, delimiter='\t')
            writer.writerows(genome_master_array)
        print('Finished with ' + genome)


master_dict_generation(pwd + '/greens')
master_df_column_gen(pwd)
# genome_tsv_generation()
        

master_tsv_generation()

