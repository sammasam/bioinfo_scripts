#!/usr/bin/env python2.7

import io
import sys

#------------------------------------------------------------------#
# Sam Bryson, Oregon State University, Dept. of Microbiology       #
# 19 August 2014                                                   #
#------------------------------------------------------------------#

# kraken_taxonomy.py

# assigns taxonomy tree to kraken output

# input three files:
    
    # -k --> kraken.output
    
    # -n --> NCBI taxonomy nodes.dmp file
    
    # -m --> NCBI taxonomy names.dmp file
    
#return outputfile:
    
    # -o --> contigs taxonomy
    # list of contigs with NCBI taxonomy
        #<contig><taxID><domain><phyla><class><sub_class><order><sub_order>...<genus><species>
    
#------------------------------------------------------------------#
# functions
#------------------------------------------------------------------#

def ParseOptions (argList):
    count = 0;
    for i in range(1, len(argList)):
        if argList[i] == '-k':
            kFile = argList[i+1]
            count += 1
        elif argList[i] == '-n':
            nodes = argList[i+1]
            count += 1
        elif argList[i] == '-m':
            names = argList[i+1]
            count += 1;
        elif argList[i] == '-o':
            outFile = argList[i+1]
            count += 1
    if count < 4:
        print("usage for: assign_kraken_taxonomy.py\n" +
              "\t-k --> kraken output file\n" +
              "\t-n --> NCBI taxonomy nodes.dmp file\n" +
              "\t-m --> NCBI taxonomy names.dmp file\n" +
              "\t-o --> output file ... \n")
        quit()
    return(kFile,nodes,names,outFile)

def LoadTaxonomy (nodes, names): #returns dict: taxid => (parent, rank, scientific name)
    print("\n<>creating dictionary for NCBI taxonomy")
    taxaDict = {}
    ranks = []
    print("\t<>loading nodes.dmp file")
    nIN = io.open(nodes)
    for line in nIN:
        lineList = line.strip().split("\t|\t")
        taxaID = lineList[0]
        parentID = lineList[1]
        rank = lineList[2]
        taxaDict[taxaID] = [parentID,rank]
        if rank not in ranks:
            ranks.append(rank)
    nIN.close()
    print("\t<>loading names.dmp file")
    mIN = io.open(names)
    for line in mIN:
        listLine = line.strip().split("\t|\t")
        taxaID = listLine[0]
        name = listLine[1]
        taxonomy = taxaDict[taxaID]
        taxonomy.append(name)
        taxaDict[taxaID] = taxonomy
        #print(str(taxaID)+"\t"+str(parentID)+"\t"+rank+"\t"+name)
    nIN.close()
    print(ranks)
    return(taxaDict)

def AssignTaxonomy (kFile, taxaDict):
    print("\n<>creating dictionary for contig taxonomy")
    kDict = {}
    print("\t<>loading kraken output file: "+kFile)
    kIN = io.open(kFile)
    for line in kIN:
        lineList = line.strip().split("\t")
        status = lineList[0]
        contig = lineList[1]
        taxaID = lineList[2]
        if status == 'C':
            taxonomy = taxaDict[taxaID]
            parentID = taxonomy[0]
            rank = taxonomy[1]
            name = taxonomy[2]
        else:
            [parentID,rank,name] = ['U','U','U']
        #print(contig+"\t"+taxaID+"\t"+rank+"\t"+name)
        toplevel = rank+"="+name
        kDict[contig] = (status,taxaID,parentID,[toplevel])
        #print(status+"\t"+taxaID+"\t"+parentID+"\t"+toplevel)
    kIN.close()
    return(kDict)

def GetLineage (kDict, taxaDict):
    cDict = {}
    for contig in kDict:
        (status,taxaID,parentID,lineage) = kDict[contig]
        if status == 'U':
            kDict[contig] = ('U', 'U')
        else:
            nodes = AddNodes(taxaDict, parentID, lineage)
            cDict[contig] = (status, nodes)
            #print(lineage)
    return(cDict)
        
def AddNodes (taxaDict, parentID, lineage):
    if parentID == '1':
        return(lineage)
    else:
        taxonomy = taxaDict[parentID]
        parentID = taxonomy[0]
        rank = taxonomy[1]
        name = taxonomy[2]
        #[parentID,rank,name] = taxaDict[parentID]
        
        currentLevel = rank+"="+name
        lineage.append(currentLevel)
        AddNodes(taxaDict, parentID, lineage)
    return(lineage)

def WriteContigTaxonomy (cDict,outFile):
    print("\t<>writing contig taxonomy to: "+outFile)
    fOUT = io.open(outFile, 'a')
    fOUT.write(unicode("#Contig\tStatus\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n"))
    for contig in cDict:
        (status, lineage) = cDict[contig]
        (domain,phylum,classs,order,family,genus,species) = ParseLineage(lineage)
        fOUT.write(unicode(contig+"\t"+status+"\t"+domain+"\t"+phylum+"\t"+
                           classs+"\t"+order+"\t"+family+"\t"+genus+"\t"+species+"\n"))
        #print(contig+"\t"+status+"\t"+domain+"\t"+phylum+"\t"+classs+"\t"+order+"\t"+family+"\t"+genus+"\t"+species)
    fOUT.close()
    return()

def ParseLineage (lineage):
    (domain,phylum,classs,order,family,genus,species) = ('U','U','U','U','U','U','U')
    if lineage == 'U':
        return(domain,phylum,classs,order,family,genus,species)
    else:
        for l in lineage:
            rank = l.split('=')
            if rank[0] == 'superkingdom':
                domain = rank[1]
            elif rank[0] == 'phylum':
                phylum = rank[1]
            elif rank[0] == 'class':
                classs = rank[1]
            elif rank[0] == 'order':
                order = rank[1]
            elif rank[0] == 'family':
                family = rank[1]
            elif rank[0] == 'genus':
                genus = rank[1]
            elif rank[0] == 'species':
                species = rank[1]
    return(domain,phylum,classs,order,family,genus,species)
        
#------------------------------------------------------------------#
# Main
#------------------------------------------------------------------#

if __name__ == '__main__':
    (kFile,nodes,names,outFile) = ParseOptions(sys.argv)
    taxaDict = LoadTaxonomy(nodes, names)
    kDict = AssignTaxonomy(kFile, taxaDict)
    cDict = GetLineage (kDict, taxaDict)
    WriteContigTaxonomy(cDict,outFile)
    
    