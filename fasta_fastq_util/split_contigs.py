#!/usr/bin/env python2.7

#------------------------------------------------------------------#
# Sam Bryson, Oregon State University, Dept. of Microbiology       #
# 5 November 2014                                                 #
#------------------------------------------------------------------#

import io
import sys

### take a fasta file input -i <fasta_input>
    # split contigs to even lengths of -l
    # output fast -f <output>

#------------------------------------------------------------------#
# functions
#------------------------------------------------------------------#


def ParseOptions (aList):
    count = 0
    for i in range(1, len(aList)):
        if aList[i] == '-i':
            inFile = aList[i+1]
            count += 1
        elif aList[i] == '-o':
            outFile = aList[i+1]
            count += 1
        elif aList[i] == '-l':
            length = int(aList[i+1])
            count += 1
    if count < 3:
        print("split contigs by length\n\tusage:\n\t\t-i <fasta_input>\n\t\t-o <fasta_output>"+
              "\n\t\t-l <length>\n")
        quit()
    return(inFile,outFile,length)

def ContigSplitter (inFile,outFile,length):
    print("splitting contigs in fasta file: " + inFile+"\n\toutput lengths of: "+str(length))
    fIN = io.open(inFile)
    fOUT = io.open(outFile, 'a')
    seq = ""
    header = ""
    for line in fIN:
        line = line.strip()
        if line[0] == ">":
            if seq:
                parts = GetParts(seq,header,length)
                fOUT.write(unicode(parts))
                seq = ""
            header = line[1:]
        else:
            seq += line
    if seq:
        parts = GetParts(seq,header,length)
        fOUT.write(unicode(parts))
        seq = ""
    fIN.close()
    fOUT.close()
    return()

def GetParts (seq,header,length):
    outString = ""
    n = len(seq)
    a = n/length
    b = (a + 1) * length
    c = b - n
    d = c/a
    ints = [d]*(a-1)
    ints.append((c-sum(ints)))
    #print(header+"\nlength = "+str(n)+"\nparts = "+str(a+1))
    #print(str(ints))
    start = 0
    end = length
    for i in range(0,(len(ints))):
        part = seq[start:end]
        outString = outString+">"+header+"."+str(i)+"\n"+part+"\n"
        start = start + length - ints[i]
        end = start + length
    part = seq[start:end]
    outString = outString+">"+header+"."+str(i+1)+"\n"+part+"\n"
    return(outString)

#------------------------------------------------------------------#
# Main
#------------------------------------------------------------------#

if __name__ == '__main__':
    (inFile,outFile,length) = ParseOptions(sys.argv)
    ContigSplitter(inFile,outFile,length)
