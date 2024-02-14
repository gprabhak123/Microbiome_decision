import argparse
import pandas as pd
import csv
import numpy as np
import os
import re
import pickle as pic
from fusion import isFusion, geneFusions, genDict
from parseXML import parse
from pprint import pprint
import os
import subprocess

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

parser = argparse.ArgumentParser(description='genome.')
parser.add_argument('--genome', type=str, default='GCF_009648975.1', help='Genome ID')
args = parser.parse_args()
if not os.path.exists(args.genome):
    print(f"Genome folder not found: {args.genome}")
    exit()

TCDB_seqs =args.genome + "/tcdb_seqs"
if os.path.exists(TCDB_seqs):
    os.system(f"rm -r {TCDB_seqs}")
os.mkdir(TCDB_seqs)
data = pd.read_csv(args.genome + '/results.tsv', sep='\t')
hit_tcid_array = data['Hit_tcid'].unique()
command_1=[f"extractTCDB.pl -i {tcid} -o {args.genome}/tcdb_seqs -f fasta" for tcid in hit_tcid_array]
joined_commands = ';'.join(command_1)

execute_command(joined_commands)
command2=f"cd {args.genome};rm -f all.faa;cat tcdb_seqs/*faa >> all.faa &&hmmtop -if=all.faa -of=hmmtop.out -sf=FAS -pi=spred -is=pseudo"
execute_command(command2)