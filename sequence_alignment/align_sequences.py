#!/usr/bin/env python

"""
    usage:
        align_sequences [options] seq1.fa seq2.fa
    where the options are:
        -h,--help : print usage and quit
        -m,--match: score of a match in the alignment [2]
        -x,--mismatch: penalty for a mismatch in the alignment [1]
        -g,--gapopen: penalty for opening a new gap [4]
        -e,--gapextend: penalty for extending a gap [1]
"""

from sys import argv, stderr
from getopt import getopt, GetoptError

# a simple function to read the name and sequence from a file
# The file is expected to have just one contig/sequence. This function
# checks the assumption and complains if it is not the case.
def read_single_contig_fasta(filename):
    names = []
    sequences = []
    with open(filename, 'r') as f:
        line = f.readline()
        assert line.startswith(">")
        names.append(line.strip().split("\t"))
        sequence = ""
        for line in f:
            if line.startswith(">"):
                sequences.append(sequence)
                names.append(line.strip().split("\t"))
                sequence = ""
            else:
                for x in line.strip():
                    if x not in ["A", "C", "G", "T"]:
                        print("Unknown nucleotide {}".format(x), file=stderr)
                        exit(3)
                sequence += line.strip()

    sequences.append(sequence)
    assert len(names) == 1
    assert len(sequences) == 1
    return names[0], sequences[0]


import numpy as np

def smith_waterman(seq1, seq2, match, mismatch, gapopen, gapextend):
    max_score = 0
    M = {}
    S = np.zeros([len(seq1) + 1, len(seq2) + 1])
    
    for i in range(0,len(seq1)+1):
        for j in range(0,len(seq2)+1):
            
            M['[' + str(i) + ', ' + str(j) + ']'] = ""
            
            if i == 0 or j == 0:
                pass
            else:
                if seq1[i-1] == seq2[j-1]:
                    s = S[i-1,j-1] + match
                else:
                    s = S[i-1,j-1] - mismatch
                d = S[i-1,j-1] + s
                
                if M['[' + str(i-1) + ', ' + str(j) + ']'] != "u":
                    p = gapopen + gapextend
                else:
                    p = gapextend            
                u = S[i-1,j] - p
                
                if M['[' + str(i) + ', ' + str(j-1) + ']'] != "l":
                    p = gapopen + gapextend            
                else:
                    p = gapextend                
                l = S[i,j-1] - p
                
                S[i,j] = max(s,u,l,0)
                
                
                if s == S[i,j]:
                    M['[' + str(i) + ', ' + str(j) + ']'] = "d"
                elif u == S[i,j]:
                    M['[' + str(i) + ', ' + str(j) + ']'] = "u"
                elif l == S[i,j]:
                    M['[' + str(i) + ', ' + str(j) + ']'] = "l"
                else:
                    M['[' + str(i) + ', ' + str(j) + ']'] = ""
                                  
                    
            if S[i,j] >= max_score:
                max_score = int(S[i,j])
                max_i = i
                max_j = j

    ## trace back
    alnseq1 = ""
    alnseq2 = ""

    m = max_i #index the base of seq1 
    n = max_j #index the base of seq2
    
    i = 0 #index the base of alnseq1
    j = 0 #index the base of alnseq2
    
    
    while 1:
        m = max_i - i
        n = max_j - j
        
        if M['[' + str(m) + ', ' + str(n) + ']'] == "d":
            alnseq1 = seq1[m-1] + alnseq1
            alnseq2 = seq2[n-1] + alnseq2
            i = i+1
            j = j+1
        if M['[' + str(m) + ', ' + str(n) + ']'] == "u":
            alnseq1 = seq1[m-1] + alnseq1
            alnseq2 = "-" + alnseq2
            i = i+1
        if M['[' + str(m) + ', ' + str(n) + ']'] == "l":
            alnseq1 = "-" + alnseq1
            alnseq2 = seq2[n-1] + alnseq2
            j = j+1
        
        if M['[' + str(m) + ', ' + str(n) + ']'] == "":
            break    

    return max_score, alnseq1, alnseq2
    

def main(filename1, filename2, match, mismatch, gapopen, gapextend):
    # read the name and sequence from the file
    name1, seq1 = read_single_contig_fasta(filename1)
    name2, seq2 = read_single_contig_fasta(filename2)

    # this function takes as input two nucleotide sequences along with
    # scores for an alignment match, mismatch, opening a new gap, and 
    # extending an existing gap. This should return the maximum alignment
    # score as well as the alignment. For examples see the testdriver script
    max_score, alnseq1, alnseq2 = smith_waterman(seq1, seq2, 
                                  match, mismatch, gapopen, gapextend)
    
    print("Maximum alignment score: {}".format(max_score))
    print("Sequence1 : {}".format(alnseq1))
    print("Sequence2 : {}".format(alnseq2))

if __name__ == "__main__":
    try:
        opts, args = getopt(argv[1:],
                     "hm:x:g:e:",
                     ["help", "match=", "mismatch=", "gapopen=", "gapextend="])
    except GetoptError as err:
        print(err)
        print(__doc__, file=stderr)
        exit(1) 

    match = 2
    mismatch = 1
    gapopen = 4
    gapextend = 1

    for o, a in opts:
        if o in ("-h", "--help"):
            print(__doc__, file=stderr)
            exit()
        elif o in ("-m", "--match"):
            match = float(a)
        elif o in ("-x", "--mismatch"):
            mismatch = float(a)
        elif o in ("-g", "--gapopen"):
            gapopen = float(a)
        elif o in ("-e", "--gapextend"):
            gapextend = float(a)
        else:
            assert False, "unhandled option"

    if len(args) != 2:
        print(__doc__, file=stderr)
        exit(2)

    main(args[0], args[1], match, mismatch, gapopen, gapextend)

