#! /usr/bin/env python

import sys
import requests

'''
A script to parse the TCDB substrate data page and return it as a dictionary.
'''

def open_page(url):

    return requests.get(url).text

def parse_substrate(substrate):

    chebi,name = substrate.split(';')
    return (chebi,name)

def get_substrate_data(url):

    substrates = {}

    contents = open_page(url).split('\n')

    for line in contents:

        if line != '':

            tcid,substrate_info = line.split('\t')

            substrates[tcid] = []

            if '|' in substrate_info:

                substrate_info = substrate_info.split('|')

                for substrate in substrate_info:

                    data = parse_substrate(substrate)
                    substrates[tcid].append(data)

            else:

                data = parse_substrate(substrate_info)
                substrates[tcid].append(data)


    return substrates



if __name__ == "__main__":

    url = 'https://tcdb.org/cgi-bin/substrates/getSubstrates.py'

    get_substrate_data(url)
