import subprocess
import os
import copy
from subprocess import call

#Perform mmseqs_createdb (convert Fasta files into mmseqs's dbtype/database structure)
def runCreateDB(fastaFileName):
    createDB_cmd = ['mmseqs', 'createdb', fastaFileName, fastaFileName + 'DB']
    if os.path.isfile(fastaFileName + 'DB') == False:
        subprocess.run(createDB_cmd, check=True)


#Perform mmseqs_search (calculate alignment from query and target databases)
def runSearch(query, target, result):
    search_cmd = ['mmseqs', 'search', query, target, result, 'tmp', '-a', '--gap-open', 'aa:9,nucl:5', '--gap-extend', 'aa:1,nucl:2']
    if os.path.isfile(result + ".dbtype") == False:
        subprocess.run(search_cmd, check = True)



#Perform mmseqs_convertalis (convert alignment result from 'Search' module to BLAST tab)
def runConvertalis(query, target, result):
    column_labels = 'query,target,alnlen,evalue,bits,pident,qstart,qend,qlen,qcov,tstart,tend,tlen,tcov,qaln,taln,qseq,tseq'
    convertalis_cmd = ['mmseqs', 'convertalis', query, target, result, result + ".tab", "--format-mode", "4","--format-output", column_labels]
    if os.path.isfile(result + ".tab") == False:
        subprocess.run(convertalis_cmd, check = True)



#Convert BLAST tab file generated by mmseqs Convertalis into python dictionary object
def blastDict(result):
    #Convert BLAST tab file into dictionary
    #Open the tab file for reading
    with open(result + ".tab", 'r') as f:
         # Read the header line and split it by tabs
        header = f.readline().strip().split('\t')
        # Initialize an empty dictionary
        mmseqs_dict = {}
        # Read through each line of the file
        for line in f:
            # Split the line by tabs
            fields = line.strip().split('\t')
            fields[2] = int(fields[2]) #alnlen
            fields[3] = float(fields[3]) #evalue
            fields[4] = int(fields[4]) #bits
            fields[5] = float(fields[5]) #pident
            fields[6] = int(fields[6]) #qstart
            fields[7] = int(fields[7]) #qend
            fields[8] = int(fields[8]) #qlen
            fields[9] = float(fields[9]) #qcov
            fields[10] = int(fields[10]) #tstart
            fields[11] = int(fields[11]) #tend
            fields[12] = int(fields[12]) #tlen
            fields[13] = float(fields[13]) #tcov

            # Create a dictionary using the header as the keys and the fields as the values
            blast_result = {header[i]: fields[i] for i in range(len(header))}
            # Use the first field as the key to the nested dictionary
            key = blast_result.pop(header[0])
            # Add the blast result to the nested dictionary
            mmseqs_dict[key] = blast_result
        
        mmseqs_dict = dict(sorted(mmseqs_dict.items()))
        return mmseqs_dict


#Create dictionary of just targets to remove duplicate target accessions
def target_OnlyDict(mmseqsDict):
    target_only = {}
    for key in mmseqsDict:
        target_only[mmseqsDict[key]['target']] = mmseqsDict[key]['tseq']
    target_only = dict(sorted(target_only.items()))
    return target_only


#Create fasta file used for Hmmtop's TMS calculation from python dictionaries (query and target's accessions and aligned sequences)
def fastaConversion(result, mmseqsDict, targetDict):
    #Convert Result From mmseqs alignment to fasta format
    if os.path.isfile(result + '.faa') == False: #check if file already exists
        call('touch ' + result+ '.faa', shell = True)

    ofile = open(result + '.faa', "w") #open file for writing 
    faaLnCount = 0

    for key in mmseqsDict: #write in query accessions and sequences 
        faaLnCount += 1 #keep track of number of iterations in mmseqs_dict in creation of the .faa file
        ofile.write(">" + key + "\n" + mmseqsDict[key]['qseq'] + "\n")

    tfaaLnCount = 0
    for key in targetDict: #write in target accesssions and sequences
        tfaaLnCount += 1 
        ofile.write(">" + key + "\n" + targetDict[key] + "\n")

    ofile.close()

    print("number of sequences contained the .faa file is:" + str(faaLnCount + tfaaLnCount)) #2020 -> 1010 sequences
    print("Converted result file from mmseqs search to FASTA format as: " + result + ".faa" + " for use in hmmtop")

    faaSizeCmd = ['wc', '-l', result + '.faa']
    subprocess.run(faaSizeCmd, check = True) #2006 -> 1003 sequences


#Run Hmmtop on a fasta file (storing accessions and aligned sequences)
def runHmmtop(fastaFile, hmmResult):
    #Perform hmmtop
    hmmtop_cmd = ['hmmtop', '-if=' + fastaFile, '-of=' + hmmResult +'.hmmtop', '-sf=FAS', '-pi=spred', '-is=pseudo']
    if os.path.isfile(hmmResult + '.hmmtop') == False:
        subprocess.run(hmmtop_cmd)


#Command that pulll substrate data from tcdb website
def getSubstrateData(fileName):
    substrate_cmd = ['wget', '-O', fileName + '.tsv', 'https://tcdb.org/cgi-bin/substrates/getSubstrates.py']
    if os.path.isfile(fileName + '.tsv') == False:
        subprocess.run(substrate_cmd, check = True)


#Record results of Hmmtop in python dictionary (accession as key; protein length, number of TMS, and list of TMS as values)
def hmmtopDict(hmmResult):
    # Open the .hmmtop file for reading
    with open( hmmResult + '.hmmtop', 'r') as f:
        # Initialize an empty dictionary
        hmmtop_dict = {}
        # Manually create names for fields
        fieldNames = ["protlen", "ntms","tms"]
        # Read through each line of the file
        for l in f:
            # Split the line by spaces
            fields = l.strip().split()
            # empty list for individual line
            hmmtop_result = {}
            # create list to store pairs of tms
            tms_list = []
            tempTMSPos = 0
            for i in range (len(fields)):
                if fields[i] == 'IN' or fields[i] == 'OUT':
                    tempTMSPos = i
            for i in range(tempTMSPos + 2, len(fields), 2):
                tms_pairs = [int(fields[i]), int(fields[i+1])]
                tms_list.append(tms_pairs)
            hmmtop_result = {fieldNames[0]:int(fields[1]), fieldNames[1]:int(fields[4]), fieldNames[2]:tms_list}
            hmmtop_dict[fields[2]] = hmmtop_result

        return hmmtop_dict


#Record substrates data in tsv file in python dictionary
def substrateDict(fileName):
    with open(fileName + '.tsv', 'r') as f:
        # Initialize an empty dictionary
        substrate_dict = {}
        # Read through each line of the file
        for line in f:
            # Split the line by tabs
            fields = line.strip().split('\t')
            chebi = fields[1].split('|')
            element = []
            for i in chebi:
                rawPair = i.split(';')
                pair = [int(rawPair[0].split(":")[1]), rawPair[1]]
                element.append(pair)
            substrate_dict[fields[0]] = element
        
        return substrate_dict


#Calculate number of TMS overlaps in each pair of aligned sequences
def overlapDict(mmseqsDict, hmmtopDict, minRes):
    overlap_dict = {}


    def oneToZeroIndexed(coords):
        newcoords = copy.deepcopy(coords)
        for i in newcoords:
            i[0] = i[0] - 1
            i[1] = i[1] - 1
        return newcoords


    #Assume coords have been already corrected to 0-based indices
    def realRMScoords (coords, alnStart, alnEnd):
        newcoords = copy.deepcopy(coords)
	    #Remove TMSs that are not in the alignment region
	    #using the equation:
	    #   relTMCScoord = hmmtopCoord - alnStart
        returnTMS = []
        for i in newcoords:
            i[0] = i[0] - alnStart
            i[1] = i[1] - alnStart

            if i [1] <= 0:
                continue

            elif i[0] >= alnEnd - alnStart:
                continue
        
            else:
                returnTMS.append(i)

        return returnTMS


    #Helper function that returns true if given starting index is within the range of a TMS, otherwise returns false
    def findTMS(givenIndex, TMSList):
        for i in TMSList:
            if i[0] <= givenIndex and i[1] >= givenIndex:
                return True
        return False
    

    #Helper function that remove a TMS from the list if given Index is within the range of a TMS
    def removeTMS(givenIndex, TMSList):
        for i in TMSList:
            if i[0] <= givenIndex and i[1] >= givenIndex:
                TMSList.remove(i) #CHECK SYNTAX
                return
        return

    
    def overlapScore (qTMSWithin, tTMSWithin, qaln, taln):
        overlaps = 0
    
        qTracking = 0
        tTracking = 0

        #copy two lists of TMSWithin, SYNTAX MIGHT NOT BE CORRECT!!!!
        tempQWithin = qTMSWithin
        tempTWithin = tTMSWithin

        residueCount = 0
        for i in range (0, len(qaln)): #aligned region has the same length for query and target (counting the gaps)

            #if encounter gaps, actual indices ("Tracking") don't increment, and no need to compare TMS
            if qaln[i] == "-":
                qTracking = qTracking - 1
            if taln[i] == "-":
                tTracking = tTracking - 1

            #identical Residue Position Found
            if findTMS(qTracking, tempQWithin) == True and qaln[i] != "-" and taln[i] != "-":
                if findTMS(tTracking, tempTWithin) == True :
                    residueCount = residueCount + 1
                else:
                    residueCount = 0 #making sure next residue count doesn't inherit previous count
            elif findTMS(qTracking, tempQWithin) == False :
                residueCount = 0 #making sure next residue count doesn't inherit previous count

            if residueCount == minRes:
                overlaps = overlaps + 1
                residueCount = 0 #making sure next residue count doesn't inherit previous count

                #If one pair TMS is already counted as overlapped, they can't be used for future comparison
                removeTMS(qTracking, tempQWithin) 
                removeTMS(tTracking, tempTWithin)

            #tracking actual indices
            qTracking = qTracking + 1
            tTracking = tTracking + 1

        #print(tempQWithin)
        #print(tempTWithin)
        return overlaps


    for key in mmseqsDict:
        qaln = mmseqsDict[key]['qaln']
        taln = mmseqsDict[key]['taln']
        tAccession = mmseqsDict[key]['target']

        # - 1 to convert from 1-indexing to 0-indexing
        qalnStartPos = mmseqsDict[key]['qstart'] - 1
        qalnEndPos = mmseqsDict[key]['qend'] - 1
        talnStartPos = mmseqsDict[key]['tstart'] - 1
        talnEndPos = mmseqsDict[key]['tend'] - 1

        rawqTMS = hmmtopDict[key]['tms']
        rawtTMS = hmmtopDict[tAccession]['tms']

        qTMS = oneToZeroIndexed(rawqTMS)
        tTMS = oneToZeroIndexed(rawtTMS)

        qTMSWithin = realRMScoords(qTMS, qalnStartPos, qalnEndPos)
        tTMSWithin = realRMScoords(tTMS, talnStartPos, talnEndPos)
        
        
        overlaps = overlapScore(qTMSWithin, tTMSWithin, qaln, taln)
        targetAndOverlap = {}
        targetAndOverlap['target'] = tAccession
        targetAndOverlap['alignedTMS'] = overlaps

        overlap_dict[key] = targetAndOverlap	
   
        
    return overlap_dict



'''
#For Testing
tempMmseqs = {}
tempData = {}
tempData['qaln'] = 'MGFDIGGDIGKPLKDAFDKFGADIKMTFLTVLNWMK--WISIG------ILIVISVI-------LICKIIKVLFQCGKCLLSCFGFCKK'
tempData['taln'] = 'MGFSINFD---PIINKFREFQTNINHNINEQLDKLKMVWINLGSHIKYWFIIIISILTILFILFLLIKITKLILNCKKIFSCCCNVCCK'
tempData['target'] = '1.A.100.1.1-B2X7D9'
tempData['qstart'] = 1
tempData['qend'] = 74
tempData['tstart'] = 22
tempData['tend'] = 107
tempMmseqs['1.A.95.2.5-YP_009361958'] = tempData

tempHmmtop = {}
tempQhmmData = {}
tempQhmmData['tms'] = [[26, 49], [54,72]]
tempThmmData = {}
tempThmmData['tms'] = [[68, 92]]
tempHmmtop['1.A.95.2.5-YP_009361958'] = tempQhmmData
tempHmmtop['1.A.100.1.1-B2X7D9'] = tempThmmData
overlap_dict = overlapDict(tempMmseqs, tempHmmtop, minRes)
'''



def main():
    qFasta_str = input("Enter query fasta file you wish to convert to mmseqs-database type:")
    tFasta_str = input("Enter target fasta file you wish to convert to mmseqs-database type:")
    runCreateDB(qFasta_str)
    runCreateDB(tFasta_str)
    
    #Ask for names of query and target databases, and Ask to name the result database
    query_str = input("Enter query database: ")
    target_str = input("Enter target database: ")
    result_str = input("Name the result database from mmseqs search: ")

    runSearch(query_str, target_str, result_str)
    runConvertalis(query_str, target_str, result_str)

    mmseqs_dict = blastDict(result_str)
    # Print the resulting dictionary for checking

    ################################
    '''
    eValueRangeDict = {}
    for key in  mmseqs_dict:
        #if mmseqs_dict[key]['evalue'] >= 1e-15 and mmseqs_dict[key]['evalue'] <= 1e-10:
        #if mmseqs_dict[key]['evalue'] > 1e-10 and mmseqs_dict[key]['evalue'] <= 1e-7:
        if mmseqs_dict[key]['evalue'] > 1e-7 and mmseqs_dict[key]['evalue'] <= 1e-4:
            eValueRangeDict[key] = mmseqs_dict[key]
    
    eValueRangeDict = dict(sorted(eValueRangeDict.items()))
    print(len(eValueRangeDict))
    '''
    ################################

    #print(mmseqs_dict)
    print("length of mmseqs_dict is:" + str(len(mmseqs_dict)))

    target_only_dict = target_OnlyDict(mmseqs_dict)
    print("length of target dictionary is:" + str(len(target_only_dict)))

    fastaConversion(result_str, mmseqs_dict, target_only_dict)

    faa_str = input('Enter file name to run hmmtop on: ')
    hmmResult_str = input('Name the result file from hmmtop (name only, without .hmmtop): ')

    runHmmtop(faa_str, hmmResult_str)

    #Download Substrate Data
    substrate_str = input('Name the file to store substrate data (without .tsv): ')
    getSubstrateData(substrate_str)

    hmmtop_dict = hmmtopDict(hmmResult_str)

    # Print the resulting dictionary
    print(hmmtop_dict)
    ##########################################
    #for key in eValueRangeDict:
        #print(hmmtop_dict[key])
    ##########################################

    print("length of the hmmtop_dict is:" + str(len(hmmtop_dict))) #1002

    substrate_dict = substrateDict(substrate_str)

    #print substrate dictionary
    #print(substrate_dict)
    #print(len(substrate_dict))

    print ("-- Calculating TMS overlap based on dictionaries generated by previously-ran mmseqs search and hmmtop --")


    minRes = 8 #minimum number of identical residues for a pair of TMS to be considered overlapping
    yn = input('Do you wish to change the minimum number of identical residues for TMS to be considered overlapping? -- default: 8 ; type y or n:  ')
    if (yn == 'y'):
        minRes = input('Enter the desired minimum residues number:')

    overlap_dict = overlapDict(mmseqs_dict, hmmtop_dict, minRes)
    #print(overlap_dict)
    #print(len(overlap_dict))

    #######################################
    '''
    Muelsyse = 0
    for key in eValueRangeDict:
        #print(key + ' : ' + str(overlap_dict[key]))
        Muelsyse = Muelsyse + 1

        
        #swTest_cmd = ['bash', 'swTest.ssh', key, eValueRangeDict[key]['target']]
        
        print(str(Muelsyse) + '. ' + key + ' : ' + str(overlap_dict[key]))
        #subprocess.run(swTest_cmd, check = True)
        
    print(Muelsyse)
    '''
    #######################################
    
# main()