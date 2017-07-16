#!/usr/bin/env python2.7

#------------------------------------------------------------------#
# gi_to_taxid.py
#
# looks up taxid #'s for gi #'s found in a table
# such as column 2 of a BLAST output
#
# Requires NCBI nodes.dmp and names.dmp files
#
# Ouput used in step 2 of the assign contig taxonomy workflow
#
# Sam Bryson, Oregon State University, Dept. of Microbiology       #
# 15 February 2016                                                 #
#------------------------------------------------------------------#

import io
import sys

# input four files:
    
    # -f --> list or BLAST tabular output with ncbi id headers ex: gi|#########|ref|YZ_######.#|
    
    # -g --> NCBI gi to taxid file: protein = gi_taxid_prot.dmp ; nucleotide = gi_taxid_nucl.dmp
    
    # -n --> NCBI taxonomy nodes.dmp file
    
#indicate which column the ncbi header is in: 0 = first column

    # -c --> column number

#return outputfile:
    
    # -o --> output_file ; list of unique gi numbers with NCBI taxid numbers D|P|C|O|F|G|S
        #<gi><d_taxID><p_taxID><c_taxID><o_taxID><f_taxID><g_taxID><s_taxID>
    
#------------------------------------------------------------------#
# functions
#------------------------------------------------------------------#

def parse_options (argList):
    count = 0
    for i in range(1, len(argList)):
        if argList[i] == '-f':
            inFile = argList[i+1]
            count += 1
        elif argList[i] == '-c':
            column = int(argList[i+1])
            count += 1
        elif argList[i] == '-t':
            gi_taxid_tree = argList[i+1]
            count += 2
            gi_taxid = 'false'
            nodes = 'false'
        elif argList[i] == '-g':
            gi_taxid = argList[i+1]
            count += 1
            gi_taxid_tree = 'false'
        elif argList[i] == '-n':
            nodes = argList[i+1]
            count += 1
            gi_taxid_tree = 'false'
        elif argList[i] == '-o':
            outFile = argList[i+1]
            count += 1
    if count < 5:
        print("usage for: assign_kraken_taxonomy.py\n" +
              "\t-f --> input file\n" +
              "\t-c --> sequence header column in tab delimited input file 0,1,...,etc\n" +
              "\t\tcolumn 1 in blast output format 6\n" +
              "\t-t --> gi_taxid tree file ... DATABASE SPECIFIC FILE\n" +
              "\t\t only use if not using -g and -n flags!!!\n" +
              "\t-g --> gi_taxid file ... make sure prot or nucl file type is correct!\n " +
              "\t-n --> NCBI taxonomy nodes.dmp file\n" +
              "\t-o --> output file ... \n")
        quit()
    if gi_taxid_tree != 'false':
        print('\nusing custom taxid tree for gi to taxid lookup')
    return(inFile,column,gi_taxid_tree,gi_taxid,nodes,outFile)

def load_gis (inFile, column):
    giDict = {}
    print("loading unique gi numbers from input file: "+inFile)
    inF = io.open(inFile)
    for line in inF:
        lineList = line.strip().split("\t")
        header = lineList[column]
        header_parts = header.split("|")
        gi = int(header_parts[1])
        if gi not in giDict:
            giDict[gi] = 'U'
    inF.close()
    gi_count = len(giDict)
    print("unique gi numbers: "+str(gi_count))
    return(giDict)

def load_database_taxonomy_tree (giDict, gi_taxid_tree):
    print("loading gi taxonomy tree: "+gi_taxid_tree)
    inF = io.open(gi_taxid_tree)
    for line in inF:
        lineList = line.strip().split("\t")
        if lineList[0] != 'GI':
            gi = int(lineList[0])
            lineage = lineList[1:]
            #lineage = [int(x) for x in lineage]
            lineage.reverse()
            if gi in giDict:
                giDict[gi] = lineage
    inF.close()
    for gi in giDict:
        if giDict[gi] == 'U':
            lineage = ['U']*7
            giDict[gi] = lineage
    return()
    
def load_gi_taxids (giDict, gi_taxid): #load taxIDs in reverse order [species-->,genus-->,...-->domain]
    print("loading gi_taxid file: "+gi_taxid)
    classified = 0
    giF = io.open(gi_taxid)
    for line in giF:
        lineList = line.strip().split("\t")
        gi = int(lineList[0])
        taxid = int(lineList[1])
        if gi in giDict:
            giDict[gi] = taxid
            classified += 1
    print("gis with taxIDs: "+str(classified))
    giF.close()
    return()
    
def load_taxonomy_tree (nodes): #returns dict: taxid => (parent_taxID, rank)
    print("creating dictionary for NCBI taxonomy tree from :"+nodes)
    taxonomyDict = {}
    taxonomyDict['U'] = ('U','U')
    nodesF = io.open(nodes)
    for line in nodesF:
        lineList = line.strip().split("\t|\t")
        taxid = int(lineList[0])
        parent_taxid = int(lineList[1])
        rank = str(lineList[2])
        taxonomyDict[taxid] = (parent_taxid,rank)
    nodesF.close()
    return(taxonomyDict)

def assign_taxonomy (giDict, taxonomyDict):
    print("building gi taxonomy tree")
    inF = io.open(inFile)
    for gi in giDict:
        taxid = giDict[gi]
        lineage = ['U']*7
        if taxid != 'U':
            if taxid in taxonomyDict:
                (parent_taxid, rank) = taxonomyDict[taxid]
                lineage = update_lineage (taxid, parent_taxid,rank, lineage, taxonomyDict)
        giDict[gi] = lineage
    inF.close()
    return()

def update_lineage (taxid, parent_taxid, rank, lineage, taxonomyDict):
    rank_index = {'species':0,'genus':1,'family':2,'order':3,'class':4,'phylum':5,'superkingdom':6}
    if rank in rank_index:
        index = rank_index[rank]
        lineage[index] = taxid
    if parent_taxid == 1:
        return(lineage)
    else:
        (xid, rank) = taxonomyDict[parent_taxid]
        lineage = update_lineage (parent_taxid, xid,rank, lineage, taxonomyDict)
    return(lineage)
           
def write_taxonomy (giDict,outFile):
    print("writing  taxonomy to: "+outFile)
    outF = io.open(outFile, 'a')
    outF.write(unicode("GI\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n"))
    for gi in giDict:
        lineage = giDict[gi] #[spe,gen,fam,odr,cla,phy,dom]
        lineage.reverse()
        tree = '\t'.join([str(t) for t in lineage])
        outF.write(unicode(str(gi)+"\t"+tree+"\n"))
    outF.close()
    return()

#------------------------------------------------------------------#
# Main
#------------------------------------------------------------------#

if __name__ == '__main__':
    (inFile,column,gi_taxid_tree,gi_taxid,nodes,outFile) = parse_options(sys.argv)
    giDict = load_gis (inFile, column)
    if gi_taxid_tree != 'false':
        load_database_taxonomy_tree (giDict, gi_taxid_tree)
    else:
        load_gi_taxids (giDict, gi_taxid)
        taxonomyDict = load_taxonomy_tree (nodes)
        assign_taxonomy (giDict, taxonomyDict)
    write_taxonomy (giDict,outFile)