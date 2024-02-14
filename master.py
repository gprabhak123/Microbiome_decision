#!/usr/bin/env python

import os
import argparse
import logging
import microbiome_main
import subprocess

def main():
    pwd = '/ResearchData/Microbiome/gblast'
    assembly_dir = '/ResearchData/Microbiome/Assemblies'
    log_file = 'log.txt'

    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s: %(message)s')

    if not os.path.isdir(pwd):
        logging.info(f"Directory does not exist: {pwd}")
        return

    # Parse command-line arguments for global thresholds
    parser = argparse.ArgumentParser(description='Master script for Microbiome analysis')
    parser.add_argument('-q', '--qcov', type=float, default=0.8, help='Query coverage threshold')
    parser.add_argument('-s', '--scov', type=float, default=0.9, help='Subject coverage threshold')
    parser.add_argument('-r', '--autored', type=float, default=0.7, help='Lowest possible coverage allowed')
    parser.add_argument('-m', '--membrane', type=int, default=5, help='Membrane proteins threshold')
    parser.add_argument('-t', '--tmsdiff', type=int, default=2, help='Lowest amount of TMSs to be a Membrane Protein')
    # parser.add_argument('-f', '--feattable', type=str, help='Path to genomic feature table')
    # parser.add_argument('-b', '--gbf', type=str, help='Path to gbff file')

    args = parser.parse_args()

    # Set global thresholds
    '''
    microbiome_main.Q_COV_THRESH = args.qcov
    microbiome_main.S_COV_THRESH = args.scov
    microbiome_main.AUTO_RED = args.autored
    microbiome_main.Membraneprotein_threshold = args.membrane
    microbiome_main.Hit_TMS_Diff = args.tmsdiff
    '''
    
    # Iterate over directories in the directory
    for genome in os.listdir(pwd):
        if not genome.startswith('GCF'):
            continue
        genome_path = os.path.join(pwd, genome)
        gbff_file = os.path.join(pwd, genome + '/*.gbff.gz')
        feature_table = os.path.join(pwd, genome + '/*feature_table.txt.gz')
        # Check if it's a directory and process it
        files_to_check = ['Red', 'Green', 'Yellow']
    
        if all(os.path.exists(os.path.join(genome_path, 'analysis', color + '.tsv')) for color in files_to_check):
            logging.info("All required files (Red.tsv, Green.tsv, Yellow.tsv) exist in the analysis folder.")
        else:
            # Set the genome global variable before calling main()
            # microbiome_main.GENOME = genome_path
            logging.info(f"Processing genome: {genome_path}")
            # Call 'microbiome_main.py' main function directly
            #microbiome_main.run()
            # microbiome_main.main()
            subprocess.run(['python', 'microbiome_main.py', '--genome', genome_path, '-s', str(args.scov), '-q', str(args.qcov), '-r', str(args.autored), '-m', str(args.membrane), '-t', str(args.tmsdiff), '--gbf', gbff_file, '--feattable', featuretable])
            print('finished with: ' + genome)
main()
