import numpy as np
import pandas as pd

df = pd.read_csv('family_matrix.tsv')

rows_with_only_one_1 = df[(df == 1).sum(axis=1) == 1]


all_1s = df[(df == 1).sum(axis=1) == 29]

rows_with_only_one_1.to_csv('Fams_found_only_once.tsv')
all_1s.to_csv('Fams_in_all_genomes.tsv')