#!/usr/bin/env python2.7

#------------------------------------------------------------------#
# assign_contig_taxid.py
#
# assigns taxonomy to contigs based on homology scores from a
# BLAST search against a reference database
#
# Requires a taxa ID file for all hits in column #2 of the
# BLAST output table
#
# Sam Bryson, Oregon State University, Dept. of Microbiology       #
# 15 February 2016                                                 #
#------------------------------------------------------------------#

import io
import sys
import numpy as np

#------------------------------------------------------------------#
# objects
#------------------------------------------------------------------#

class Contig():
    def __init__(self, contig_name):
        self.contig_name = contig_name
        #self.contig_sequence = contig_sequence
        #self.contig_length = 0
        self.contig_cds = []
        self.contig_cds_count = 0
        self.contig_taxa = {}
        self.contig_taxonomy = {}
        self.ranks = ['domain','phylum','class','order','family','genus','species']
        for rank in self.ranks:
            self.contig_taxa[rank] = {}
            self.contig_taxonomy[rank] = []
        self.taxonomy = []
        self.scores = []
        
    def add_CDS (self, CDS):
        self.contig_cds.append(CDS)
        self.contig_cds_count = len(self.contig_cds)
        
    def get_contig_taxonomy (self):
        for rank in self.ranks:
            for cds in self.contig_cds:
                (name, score) = (cds.lca[rank], cds.score)
                if name in self.contig_taxa[rank]:
                    self.contig_taxa[rank][name] += score
                else:
                    self.contig_taxa[rank][name] = score
        
    def get_max_scores (self):
        for rank in self.ranks:
            max_score = 0.0
            best_taxa = ""
            total_score = 0.0
            for name in self.contig_taxa[rank]:
                score = self.contig_taxa[rank][name]
                total_score += score
                if score > max_score:
                    if name == 'U':
                        u_score = score
                        u_taxa = 'U'
                    else:
                        max_score = score
                        best_taxa = name
            if max_score > 0.0:
                self.taxonomy.append(best_taxa)
                if total_score > 0.0:
                    self.scores.append(max_score/total_score)
                else:
                    self.scores.append(0.0)
            else:
                self.taxonomy.append(u_taxa)
                if total_score > 0.0:
                    self.scores.append(u_score/total_score)
                else:
                    self.scores.append(0.0)
                    
                
#------------------------------------------------------------------#

class CDS():
    def __init__(self, cds_name, cds_contig):
        self.cds_name = cds_name
        self.cds_contig = cds_contig
        #self.cds_sequence = ""
        #self.cds_length = 0
        self.cds_blast_hits = []
        self.cds_taxonomy = {}   # => (dom,phy,cla,odr,fam,gen,spe)
        self.ranks = ['domain','phylum','class','order','family','genus','species']
        for rank in self.ranks:
            self.cds_taxonomy[rank] = []
        self.lca = {}
        self.score = 0.0
        
    def add_blast_hit (self, blast_hit):
        self.cds_blast_hits.append(blast_hit)
        for rank in self.ranks:
            hits = self.cds_taxonomy[rank][:]
            hits.append(blast_hit.hit_taxonomy[rank])
            self.cds_taxonomy[rank] = hits
            #print("\t"+rank+" | "+str(self.cds_taxonomy[rank]))
    
    def get_LCA_score (self, weight):
        for rank in self.ranks:
            hits = self.cds_taxonomy[rank][:]
            taxa = [h[0] for h in hits]
            taxa = set(taxa)
            self.score = np.mean([h[1] for h in hits])
            #
            if weight > 0.0:
                self.score = self.score/(weight * len(self.cds_blast_hits))
            # Divide socre by number of hits for cds > threshold score --> weight unique hits more!
            #
            if len(taxa) > 1:
                taxa = 'U'
            else:
                taxa = list(taxa)
                taxa = taxa[0]
            self.lca[rank] = taxa

#------------------------------------------------------------------#
        
class BLAST_hit():  # only create object if best hit or bitscore > x% of best hit bitscore
    def __init__(self, query_cds, hit, alignment_identity, alignemt_length, bit_score, score_type='alignment'):
        # optional score types => alignment(%identity*alignment length) or bitscore
        self.query_cds = query_cds
        self.hit = hit
        self.alignment_identity = alignment_identity
        self.alignemt_length = alignemt_length
        self.bit_score = bit_score
        self.hit_taxonomy = {}
        self.ranks = ['domain','phylum','class','order','family','genus','species']
        if score_type == 'alignment':
            self.score = self.alignment_identity * self.alignemt_length /100.0
        if score_type == 'bitscore':
            self.score = bit_score
        if score_type == '1':
            self.score = 1

    def get_hit_taxonomy(self, hit_taxonomy_dict):
        if self.hit in hit_taxonomy_dict:
            taxonomy = hit_taxonomy_dict[self.hit] # dom,phy,cla,odr,fam,gen,spe
        else:
            taxonomy = ['U','U','U','U','U','U','U']
        for r in range(len(self.ranks)):
            self.hit_taxonomy[self.ranks[r]] = (taxonomy[r], self.score)
           #print("\t"+self.query_cds+" | "+str(self.hit_taxonomy[self.ranks[r]]))

#------------------------------------------------------------------#
# functions
#------------------------------------------------------------------#

def parse_options (aList):
    count = 0
    score_type='bitscore'
    threshold = 0.9
    weight = 1.0
    min_score = ''
    for i in range(1, len(aList)):
        if aList[i] == '-i':
            blastFile = aList[i+1]
            count += 1
        elif aList[i] == '-t':
            taxaFile = aList[i+1]
            count += 1
        elif aList[i] == '-o':
            outFile = aList[i+1]
            count += 1
        elif aList[i] == '-s':
            score_type = aList[i+1]
        elif aList[i] == '-d':
            threshold = float(aList[i+1])  # Percent of best hit's bit score for inclusion of other hits in analysis
        elif aList[i] == '-w':
            weight = float(aList[i+1])  # Weight of score for CDS with multiple hits above threshold in contig assignment
        elif aList[i] == '-m':
            min_score = float(aList[i+1])  # min score to consider a hit: alignment AAI(%) or bit_score
    if count < 3:
        print("Add taxonomy data to BLAST output" +
              "\n\tusage:\n\t-i --> BLAST tabular file" +
              "\n\t-t --> taxonomy file for column #2 in BLAST file" +
              "\n\t-s --> score_type, 'alignment' (0-100 : Default) or 'bitscore' or 1 for equal weighting of CDS" +
              "\n\t-d --> threshold, percent of best hit score for secondary hit inclusion\n\t\tdefault: 0.9" +
              "\n\t-w --> weight, fraction of hit counts to use for reducing score of cds with multiple hits > threshold\n\t\tdefault:0.0"
              "\n\t-o -->output\n")
        quit()
    return(blastFile,taxaFile,score_type,threshold,weight,min_score,outFile) #ADD OPTIONS

def load_taxonomy (taxaFile):
    giDict = {}
    print("loading BLAST hits taxonomy: " + taxaFile)
    fIN = io.open(taxaFile)
    for line in fIN:
        l = line.strip().split('\t')
        l = [str(n) for n in l]
        [gi,dom,phy,cla,odr,fam,gen,spe] = l
        giDict[gi] = [dom,phy,cla,odr,fam,gen,spe]
    fIN.close()
    return(giDict)

def BLAST_taxonomy (giDict, blastFile, threshold, score_type, weight, min_score):
    print("parsing BLAST hits: " + blastFile)
    contig_annotations = {}
    current_contig = ""
    current_cds = ""
    fIN = io.open(blastFile)
    for line in fIN:
        l = line.strip().split('\t')
        pro_name = l[0]
        hit = get_gi(l[1])
        bit = float(l[11])
        aai = float(l[2])
        alen = float(l[9])
        header = pro_name.split('_')
        contig_name = '_'.join(header[:-1])
        current_hit = BLAST_hit(pro_name, hit, aai, alen, bit, score_type=score_type)
        current_hit.get_hit_taxonomy(giDict)
        if score_type == 'alignment':
            score = aai * alen /100.0
        if score_type == 'bitscore':
            score = bit
        if score_type == '1':
            score = 1
        #print(contig_name)
        #
        # if score >= min_score --> set threshold for blast hits?
        #
        if contig_name == current_contig:
            if current_cds == pro_name:
                if score >= threshold*best_hit:
                    new_cds.add_blast_hit(current_hit)
            else:
                if current_cds:
                    new_cds.get_LCA_score(weight)
                    new_contig.add_CDS(new_cds)
                current_cds = pro_name
                new_cds = CDS(pro_name,contig_name)
                new_cds.add_blast_hit(current_hit)
                best_hit = score
        else:
            if current_contig:
                if current_cds:
                    new_cds.get_LCA_score(weight)
                    new_contig.add_CDS(new_cds)
                    new_contig.get_contig_taxonomy()
                    new_contig.get_max_scores()
                    #print(new_contig.contig_name+"\tcds|"+str(new_contig.contig_cds_count)+"\n\t"+str(new_contig.taxonomy)+"\n\t"+str(new_contig.scores))
                    contig_annotations[new_contig.contig_name] = (new_contig.contig_cds_count, new_contig.taxonomy, new_contig.scores)
            current_contig = contig_name
            new_contig = Contig(contig_name)
            current_cds = pro_name
            new_cds = CDS(pro_name,contig_name)
            new_cds.add_blast_hit(current_hit)
            best_hit = score
    fIN.close()
    if current_contig:
        if current_cds:
            new_cds.get_LCA_score(weight)
            new_contig.add_CDS(new_cds)
            new_contig.get_contig_taxonomy()
            new_contig.get_max_scores()
            #print(new_contig.contig_name+"\tcds|"+str(new_contig.contig_cds_count)+"\n\t"+str(new_contig.taxonomy)+"\n\t"+str(new_contig.scores))
            contig_annotations[new_contig.contig_name] = (new_contig.contig_cds_count, new_contig.taxonomy, new_contig.scores)
    return(contig_annotations)

def get_classification_stats (contig_annotations):
    print("calculating clasification stats")
    ranks = ['Domain','Phylum','Class','Order','Family','Genus','Species']
    classified = [0,0,0,0,0,0,0]
    total_contigs = len(contig_annotations)
    print("Total contigs: "+str(total_contigs))
    for contig in contig_annotations:
        (cds_count, taxonomy, scores) = contig_annotations[contig]
        for i in reversed(range(len(taxonomy))):
            if taxonomy[i] != 'U':
                for n in range(i+1):
                    classified[n] += 1
                break
    for r in range(len(ranks)):
        rank = ranks[r]
        c = 100.0*classified[r]/total_contigs
        print("\t"+rank+": "+str(c))
    return()
    
def write_contig_taxonomy (contig_annotations, outFile):
    print("writing taxonomy to: " + outFile)
    fOUT = io.open(outFile, 'a')
    fOUT.write(unicode("ContigID\tCDS_count\tDomain; Phylum; Class; Order; Family; Genus; Species\n"))
    for contig in contig_annotations:
        (cds_count, tax, scores) = contig_annotations[contig]
        sc = [int((100*x)+0.5) for x in scores]
        fOUT.write("%s\t%d\t%s (%d); %s (%d); %s (%d); %s (%d); %s (%d); %s (%d); %s (%d)\n" % \
                   (contig,int(cds_count),tax[0],sc[0],tax[1],sc[1],tax[2],sc[2],tax[3],sc[3],\
                    tax[4],sc[4],tax[5],sc[5],tax[6],sc[6]))
    fOUT.close()
    return()

def get_gi (seq_id):
    a = seq_id.split('|')
    return(str(a[1]))

#------------------------------------------------------------------#
# main
#------------------------------------------------------------------#
    
if __name__ == '__main__':
    blastFile,taxaFile,score_type,threshold,weight,min_score,outFile = parse_options(sys.argv)
    taxaDict = load_taxonomy(taxaFile)
    contig_annotations = BLAST_taxonomy(taxaDict, blastFile, threshold, score_type, weight, min_score)
    get_classification_stats(contig_annotations)
    write_contig_taxonomy(contig_annotations, outFile)
    
    
    