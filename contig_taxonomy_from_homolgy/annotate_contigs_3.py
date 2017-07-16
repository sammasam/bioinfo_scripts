#!/usr/bin/env python2.7

#------------------------------------------------------------------#
# annotate_contigs.py
#
# assigns taxonomy names to contigs based on taxa IDs 
# determined by the assign_contig_taxid.py script
#
# Annotates up to highest taxonomic rank where the consensus
# homology score is greater than a user defined threshold score
#
# Requires NCBI nodes.dmp and names.dmp files
#
# Sam Bryson, Oregon State University, Dept. of Microbiology       #
# 15 February 2016                                                 #
#------------------------------------------------------------------#


import io
import sys

# input four files:
    
    # -c --> contig taxID file
    
    # -t --> threshold for annotation
    
    # -n --> NCBI taxonomy nodes.dmp file
    
    # -m --> NCBI taxonomy names.dmp file
    

#return outputfile:
    
    # -o --> output_file ; list of contig IDs with NCBI taxonomy
        #<contigID><taxID><domain><phyla><class><order><family><genus><species>
    
#------------------------------------------------------------------#
# functions
#------------------------------------------------------------------#

def parse_options (argList):
    count = 0
    for i in range(1, len(argList)):
        if argList[i] == '-c':
            contig_taxids = argList[i+1]
            count += 1
        elif argList[i] == '-t':
            threshold = int(argList[i+1])
            count += 1
        elif argList[i] == '-n':
            nodes = argList[i+1]
            count += 1
        elif argList[i] == '-m':
            names = argList[i+1]
            count += 1
        elif argList[i] == '-o':
            outFile = argList[i+1]
            count += 1
    if count < 5:
        print("usage for: assign_kraken_taxonomy.py\n" +
              "\t-c --> input file contig_taxids\n" +
              "\t-t --> assignment threshold percent integer (20, 50 etc)\n " +
              "\t-n --> NCBI taxonomy nodes.dmp file\n" +
              "\t-m --> NCBI taxonomy sci_names.dmp file\n" +
              "\t\tnecessary to grep <scientific name> from names.dmp --> sci_names.dmp\n" +
              "\t-o --> output file ... \n")
        quit()
    return(contig_taxids,threshold,nodes,names,outFile)

def load_taxonomy (nodes, names): #returns dict: taxid => (parent, rank, scientific name)
    print("creating dictionary for NCBI taxonomy")
    taxonomyDict = {}
    print("loading nodes.dmp file")
    nodesF = io.open(nodes)
    for line in nodesF:
        lineList = line.strip().split("\t|\t")
        taxid = int(lineList[0])
        parent_taxid = int(lineList[1])
        rank = lineList[2]
        taxonomyDict[taxid] = [parent_taxid,rank]
    nodesF.close()
    print("loading names.dmp file")
    namesF = io.open(names)
    for line in namesF:
        listLine = line.strip().split("\t|\t")
        taxid = int(listLine[0])
        name = listLine[1]
        taxonomy = taxonomyDict[taxid]
        taxonomy.append(name)
        taxonomyDict[taxid] = taxonomy
    namesF.close()
    return(taxonomyDict)

def assign_taxonomy (inFile, taxonomyDict, threshold):
    contig_taxonomy = {}
    contig_count = 0
    classified = [0,0,0,0,0,0,0]
    ranks = ['species','genus','family','order','class','phylum','domain']
    print("loading input file: "+inFile)
    print("creating dictionary for sequence taxonomy")
    inF = io.open(inFile)
    for line in inF:
        line = line.strip().split("\t")
        if line[0] != 'ContigID':
            contig_count += 1
            contig_id = line[0]
            annotation = ['U']*7
            ann_id = 'U'
            lineage = line[2]
            lineage = lineage.split('; ')
            lineage.reverse()
            for l in range(len(lineage)):
                part = lineage[l]
                [taxid, score] = part.split(' ')
                rank = ranks[l]
                score = score[1:-1]
                score = int(score)
                if score >= threshold:
                    if taxid != 'U':
                        taxid = int(taxid)
                        if taxid in taxonomyDict:
                            [parent_taxid, rank, name] = taxonomyDict[taxid]
                            annotation = update_annotation (taxid, parent_taxid,rank, name, annotation, taxonomyDict)
                            ann_id = str(taxid)
                            for n in range(l,len(classified)):
                                classified[n] += 1
                    break
            contig_taxonomy[contig_id] = (ann_id, annotation)
    inF.close()
    get_classification_stats (contig_count, classified)
    return(contig_taxonomy)

def get_classification_stats (contig_count, classified):
    print("calculating clasification stats")
    ranks = ['species','genus','family','order','class','phylum','domain']
    print("total contigs: "+str(contig_count))
    for r in range(len(ranks)):
        rank = ranks[r]
        c = 100.0*classified[r]/contig_count
        print("\t"+rank+": "+str(c))
    return()

def update_annotation (taxid, parent_taxid, rank, name, annotation, taxonomyDict):
    rank_index = {'species':0,'genus':1,'family':2,'order':3,'class':4,'phylum':5,'superkingdom':6}
    if rank in rank_index:
        index = rank_index[rank]
        annotation[index] = name
    if parent_taxid == 1:
        return(annotation)
    else:
        [xid, xrank, xname] = taxonomyDict[parent_taxid]
        annotation = update_annotation (parent_taxid, xid, xrank, xname, annotation, taxonomyDict)
    return(annotation)

def write_contig_annotations (contig_taxonomy,outFile):
    print("writing contig taxonomy annotations to: "+outFile)
    outF = io.open(outFile, 'a')
    outF.write(unicode("ContigID\tTaxID\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n"))
    for contig_id in contig_taxonomy:
        (ann_id, annotation) = contig_taxonomy[contig_id]
        annotation.reverse()
        ann_string = '\t'.join(annotation)
        outF.write(unicode(contig_id+"\t"+ann_id+"\t"+ann_string+"\n"))
    outF.close()
    return()

#------------------------------------------------------------------#
# Main
#------------------------------------------------------------------#

if __name__ == '__main__':
    (inFile,threshold,nodes,names,outFile) = parse_options (sys.argv)
    taxonomyDict = load_taxonomy (nodes, names)
    contig_taxonomy = assign_taxonomy (inFile, taxonomyDict, threshold)
    write_contig_annotations (contig_taxonomy,outFile)
    
    