#all imports of chebi_parser, comparedTransportomeAssignment, and substrate
import libchebipy
import argparse
import os,sys
import sys
import requests
import urllib
import csv
import obonet
import networkx

from decimal import Decimal
from substrate import get_substrate_data 
#from libchebipy import ChebiEntity
#from chebi_parser import find_predecessor,find_role,get_primary
from compareTransportomeAssignments import readLines


#@Parameter: substrate data, and list of tcids wish to obtain ChebiIDs from
#@Output: dictionary of TCIDs as keys and corresponding,stripped ChebiIDs (only the number part) as values
def get_chebi_id(substrateData, tcids):
    chebiIDs = {}
    for i in tcids:
        for key in substrateData:
            if key == i:
                chebiIDs[key]=[item[0] for item in substrateData[key]]
        if i not in chebiIDs:
            chebiIDs[i] = 'none'
    return chebiIDs

#@Parameter: input(dictionary with tcids as keys), ouput directly and output file name
def print_Table(inputDict, outdir, output):

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    outputFile = open('{}/{}'.format(outdir,output),'w')
    outputFile.write('#TCID\tCE\tRole\n')
    
    for key, values in inputDict.items():
        ce_list = ','.join(values['CE'])
        role_list = ','.join(values['Role'])
        
        
        outputFile.write('{}\t{}\t{}\n'.format(key, ce_list, role_list))

#Parse in the MasterTable csv file (compare group)
def parse_csv_file(file_path, delimiter=','):
    try:
        with open(file_path, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=delimiter)
            header = next(reader)  # Read the header row
            
            data = {}

            for row in reader:
                key = row[0]
                data[key] = {}
                value1 = row[2]
                value2 = row[3]
                
                if key not in data:
                    data[key] = []

                data[key]['CE'] = value1.split(",")
                data[key]['Role'] = value2.split(",")

            return data

    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

#check for potentially different entries in two dictionaries
def compare_keys(compareDict, myDict):
    diffNotIn = 0
    for key in compareDict:
        if key not in myDict:
            diffNotIn += 1
    
    return diffNotIn

#check for potentially different dictionary entries
def compare_values(compareDict, myGroup):

    entryNumDiffCE = 0
    entryNumDiffRole = 0
    entryDiffCE = 0
    entryDiffRole = 0
    minorCationLabelDiffCE = 0
    minorAnionLabelDiffCE = 0
    
    for key in compareDict:
        #print(key)

        #check for difference in the number of entries for each sub-key(CE and role)
        if len(compareDict[key]['CE']) != len(myGroup[key]['CE']):
            entryNumDiffCE +=1
            #print(compareDict[key]['CE'])
            #print(myGroup[key]['CE'])
            if len(compareDict[key]['Role']) != len(myGroup[key]['Role']):
                entryNumDiffRole +=1

        elif len(compareDict[key]['Role']) != len(myGroup[key]['Role']):
            entryNumDiffRole +=1

        else:
            #print(len(compareDict[key]['CE']))
            #print(len(compareDict[key]['Role']))
            
            for i in range(len(compareDict[key]['CE'])):
                if compareDict[key]['CE'][i].strip() != myGroup[key]['CE'][i].strip():
                    #print(compareDict[key]['CE'][i] + ' vs. ' + myGroup[key]['CE'][i])
                    
                    #cation in compare group marked as inorganic/organic cation in updated table
                    if (compareDict[key]['CE'][i].strip().split('-')[0] == 'cation(CHEBI:36916)' and myGroup[key]['CE'][i].strip().split('-')[0] == ('inorganic cation(CHEBI:36915)' or 'organic cation(CHEBI:25697)')) and (compareDict[key]['CE'][i].strip().split('-')[1] == myGroup[key]['CE'][i].strip().split('-')[1]):
                        minorCationLabelDiffCE += 1

                    elif (compareDict[key]['CE'][i].strip().split('-')[0] == 'anion(CHEBI:22563)' and myGroup[key]['CE'][i].strip().split('-')[0] == ('inorganic anion(CHEBI:24834)' or 'organic anion(CHEBI:25696)')) and (compareDict[key]['CE'][i].strip().split('-')[1] == myGroup[key]['CE'][i].strip().split('-')[1]):
                        minorAnionLabelDiffCE += 1

                    elif myGroup[key]['CE'][i].strip().split('-')[0] == 'None(None)' and compareDict[key]['CE'][i].strip().split('-')[0] == compareDict[key]['CE'][i] and compareDict[key]['CE'][i] == myGroup[key]['CE'][i].strip().split('-')[1]:
                        noneLabelCE += 1
                    
                    
                    else: 
                        print ('Sampletable: ' + str(compareDict[key]['CE'][i]))
                        print('Mytable: ' + str(myGroup[key]['CE'][i]) + '\n')

                    entryDiffCE += 1
            
            for j in range(len(compareDict[key]['Role'])):
                if compareDict[key]['Role'][j].strip() != myGroup[key]['Role'][j].strip():
                    #print(compareDict[key]['Role'][j] + ' vs. ' + myGroup[key]['Role'][j])
                    
                    '''
                    print ('1: ' + str(compareDict[key]['Role'][j]))
                    print('2: ' + str(myGroup[key]['Role'][j]))
                    '''

                    entryDiffRole += 1
            

    returnString = 'different amount of entries in CE: ' + str(entryNumDiffCE) + '\n' + 'different amount of entries in Role: ' + str(entryNumDiffRole) + '\n' + 'same amount, different content in CE: ' + str(entryDiffCE) + '\n' + 'same amount, different content in Role: ' + str(entryDiffRole) + '\n' + 'minor cation notation difference in CE: ' + str(minorCationLabelDiffCE) + '\n' + 'minor anion notation difference in CE: '+ str(minorAnionLabelDiffCE) + '\n'
    return returnString

#@Parameter: name/path to the .obo file format that needs to be parsed
'''
#@Output: graph contained by the file, dictionaries that
    1) (3 dictionaries) map primary chebi ids to their substrane name, immediate parent, and immediate relationship
    2) map substrate names to their primary chebi ids
    3) map secondary chebi ids to their primary chebi ids
'''
def oboParse(obofile):
    graph = obonet.read_obo(obofile)

    #dictionaries that maps CHEBI id(primary) to name, and name to CHEBI ID
    id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
    name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}

    #dictionary that maps primary id to secondary id(s)
    #id_to_altID = {id_: data.get('alt_id') for id_, data in graph.nodes(data=True)}

    #dictionary that maps primary id to immediate parent
    id_to_parent = {id_: data.get('is_a') for id_, data in graph.nodes(data=True)}
    
    #dictionary that maps primary id to relationship
    id_to_relationship = {id_: data.get('relationship') for id_, data in graph.nodes(data=True)}
    
    #dictionary that maps secondary id to primary id
    altID_to_id = {}
    for id_, data in graph.nodes(data=True):
        if 'alt_id' in data:
            for secondary in data['alt_id']:
                altID_to_id[secondary] = id_

    return graph,id_to_name,name_to_id,id_to_parent,altID_to_id,id_to_relationship

#Use one of the dictionaries generated by oboParse() to return the primary id to a given chebi ID
def findPrimary(secondary, altID_to_id):
    
    if secondary in altID_to_id.keys():
        return altID_to_id[secondary]
    else: #if the id given is already the primary id
        return secondary

#Use one of the dictionaries generated by oboParse() to return the substrate name of a given (primary) chebi ID
def getSubstrateName(id,id_to_name):
    return id_to_name[id]

#Use dict generated by terminalPredecessorParse() to find the desired predecessor(s) for a given chebi id
#@parameter: chebi id, ancestor dict, graph generated by oboParse()
def findPredecessor(id, ancdict, graph):
    #dfs_predecessors() return all nodes that can be traced upwards from a given node in tree
    allPredecessors = networkx.dfs_predecessors(graph,id)
    returnPredecessors = []

    #if the id itself is a terminal predecessor node
    if findPrimary(id, altIDToID) in ancdict.keys():
        returnPredecessors.append(getSubstrateName(id, idToName) + '(' + id + ')')

    #if any of all predecessors are on the terminal predecessor list
    for ancestor in allPredecessors:
        ancestor = findPrimary(ancestor, altIDToID)
        if ancestor in ancdict.keys():
            returnPredecessors.append(getSubstrateName(ancestor,idToName) + '(' + ancestor + ')')

    return returnPredecessors

#Use dict generated by oboParse() to find the desired role(s) for a given chebi id
def findRole(id, id_to_relationship, classes = None, flag = True):
    #Special case for detergent: dr. Saier wants all three roles of different hierarchy: detergent -> surfactant -> emulsifier
    if id == 'CHEBI:27780':
        return[getSubstrateName('CHEBI:35195', idToName) + '(CHEBI:35195)', getSubstrateName(id, idToName) + '(' + id + ')', getSubstrateName('CHEBI:63046', idToName) + '(CHEBI:63046)']
    
    #flag == True: immediate call to findRole
    if id == None or id_to_relationship[id] == None and flag == True:
        # in these four cases specified, the entities themselves are under the role categroy, which means their role(s) is within their predecessors
        if id != None and id == 'CHEBI:26115' or id == 'CHEBI:26672' or id == 'CHEBI:27026' or id == 'CHEBI:72316':
            id = findPrimary(id, altIDToID)
            allPredecessors = networkx.dfs_predecessors(graph,id)
            role = []
            if id in classes:
                role.append(getSubstrateName(id, idToName) + '(' + id + ')')
            for validTerminal in classes:
                if validTerminal in allPredecessors:
                    role.append(getSubstrateName(validTerminal, idToName) + '(' + validTerminal + ')')
            return role
        #id itself is a valid role to stop on (also included in the previous if clause)
        elif id != None and findPrimary(id, altIDToID) in classes:
            return [getSubstrateName(id, idToName) + '(' + id + ')']
        else: #id == None
            return []
    id = findPrimary(id, altIDToID)
    if idToParent[id] == None and flag == False:
        return[]
    
    #id itself is a valid role to stop on
    if id in classes:
        return [getSubstrateName(id, idToName) + '(' + id + ')']
    
    roles = []
    if flag == True: 
        for entry in id_to_relationship[id]:
            if entry.split(' ')[0] == 'has_role':
                targetID = findPrimary(entry.split(' ')[1], altIDToID)
                
                if targetID in classes:
                    roles.append(getSubstrateName(targetID, idToName) + '(' + targetID +')')
                
                elif targetID in id_to_relationship:
                    roles.extend(findRole(targetID, id_to_relationship, classes, False))
    #flag false -> a traversal to the relationship: has_role list has already happened
    #              -> traverse in the predecessor tree instead of role tree now
    elif flag == False:
        allPredecessors = networkx.dfs_predecessors(graph,id)
        for validTerminal in classes:
            if validTerminal in allPredecessors:
                return[getSubstrateName(validTerminal, idToName) + '(' + validTerminal + ')']
    roles = list(set(roles))
    return roles
    
#Parse the file containing desired terminal predecessor nodes, provided by Arturo (for ce ontology)
def terminalPredecessorParse(tsvFile, delimiter=' '):
    try:
        with open(tsvFile) as file:
            tsv_file = csv.reader(file, delimiter="\t")
            header = next(tsv_file)
            data = {}

            for row in tsv_file:
                category = row[0]
                entityName = row[1]
                entityIDs = []
                for id in row[2].split('/'):
                    entityIDs.append(findPrimary('CHEBI:'+ id, altIDToID))
                for primaryid in entityIDs:
                    data[primaryid] = {}
                    data[primaryid]['category'] = category
                    data[primaryid]['entityName'] = entityName
                          
            return data

    except FileNotFoundError:
        print(f"File not found: {tsvFile}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

#Record a set of tuples in a txt file
def set_of_tuples_to_txt(input_set, file_path):
    try:
        with open(file_path, 'w') as file:
            for tpl in input_set:
                line = ', '.join(map(str, tpl)) + '\n'
                file.write(line)
        print(f"Set of tuples successfully written to {file_path}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    compareDict = parse_csv_file('MasterTable.csv',',')
    substrate_dict = get_substrate_data('https://tcdb.org/cgi-bin/substrates/getSubstrates.py')
    
    #sample input (list of TCIDs)
    #input = ['1.A.1.13.3', '1.A.1.13.4', '1.A.1.13.7', '1.A.1.13.8', '1.A.1.13.10', '1.A.1.24.3', '1.A.1.28.7', '1.A.8.3.1', '1.A.11.1.4', '1.A.11.1.5']
    input = []
    for key in compareDict:
        input.append(key)
    cleanChebiIDs = get_chebi_id(substrate_dict, input)
    
    role_classes = set(['CHEBI:23888','CHEBI:33281','CHEBI:26672','CHEBI:31432','CHEBI:33229','CHEBI:23357','CHEBI:25212',
                        'CHEBI:23924', 'CHEBI:27780', 'CHEBI:35703', 'CHEBI:37958', 'CHEBI:38161', 'CHEBI:71338', 'CHEBI:62488',
                        'CHEBI:33280', 'CHEBI:25728'])

    graph, idToName, nameToId, idToParent, altIDToID, idToRelationship = oboParse('chebi.obo')

    predecessorDict = terminalPredecessorParse('substrateTypes.tsv')
    print(findPredecessor('CHEBI:17996',predecessorDict,graph))
    
    
    myGroup = {}
    noPredecessor = set()
    noRoles = set()
    for key in cleanChebiIDs:
        if cleanChebiIDs[key] == 'none':
            myGroup[key] = {}
            myGroup[key]['CE'] = ['none']
            myGroup[key]['Role'] = ['none']
        else:
            myGroup[key] = {}
            myGroup[key]['CE'] = []
            myGroup[key]['Role'] = []
            for id in cleanChebiIDs[key]:
                id = findPrimary(id, altIDToID)
                myGroup[key]['CE'].append(str(findPredecessor(id,predecessorDict,graph)) + '-' + getSubstrateName(id, idToName) + '(' + id + ')')
                if len(findPredecessor(id,predecessorDict,graph)) == 0:
                    noPredecessor.add((getSubstrateName(id,idToName), id))
                myGroup[key]['Role'].append(str(findRole(id,idToRelationship,role_classes)) + '-' + getSubstrateName(id, idToName) + '(' + id + ')')
                if len(findRole(id,idToRelationship,role_classes)) == 0:
                    noRoles.add((getSubstrateName(id,idToName), id))
                '''
                predecessor = str(findPredecessor(id, idToParent, classes=ce_classes)[0])
                predecessorID = str(findPredecessor(id,classes=ce_classes)[1])
                CE = str(chebi_entity.get_name())
                CEID = str(chebi_entity.get_id())
                role = str(find_role(i,classes=role_classes)[0])
                roleID = str(find_role(i,classes=role_classes)[1])
                #Predecessor(Predecessor Chebi ID)-CE name(CE Chebi ID)
                if predecessor == 'None' and predecessorID =='None':
                    myGroup[key]['CE'].append(CE + '(' + CEID + ')')
                else:
                    myGroup[key]['CE'].append(predecessor + '(' + predecessorID +')' +'-'+ CE + '(' + CEID + ')')
                #Role(Role Chebi ID)-CE name(CE Chebi ID)
                if role == 'None' and roleID == 'None':
                    myGroup[key]['Role'].append(CE + '(' + CEID + ')')
                else:
                    myGroup[key]['Role'].append(role + '(' + roleID +')' +'-'+ CE + '(' + CEID + ')')
                '''



    #Both compareDict and myGroup have a size of 1070 keys
   
    #print(compare_values(compareDict,myGroup))
    #Sample output directory and output file name
                
    print_Table(myGroup, '/Users/siyanlin/Downloads', 'FinalTable.tsv')
    set_of_tuples_to_txt(noRoles,'/Users/siyanlin/Downloads/NEWnoRoles3.txt')
    set_of_tuples_to_txt(noPredecessor, '/Users/siyanlin/Downloads/noPredecessors2.txt')
    #print(readLines('MasterTable.csv'))

#Finding: 1.A.1.28.7 not in output file, Likely due to it not being in substrate data

    