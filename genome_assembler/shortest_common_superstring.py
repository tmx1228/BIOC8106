#!/usr/bin/env python

"""
    usage:
        shortest_common_superstring.py [options] sequences.txt
    where the options are:
        -h,--help : print usage and quit

    sequences.txt is a file with one sequence to be included in each row.
"""

from sys import argv, stderr
from getopt import getopt, GetoptError
import itertools

def read_sequences(filename):
    '''Read the sequences (one per line) from the filename and return a list
    '''
    sequences = []
    with open(filename, 'r') as f:
        for line in f:
            sequences.append(line.strip())
    return sequences

def calculate_scs(reads):
    """
    Implement the greedy shortest-common-superstring strategy discussed in 
    class. From the reads, find two string with the maximal overlap and merge 
    them. Keep doing this until you have only 1 string left
    """
    scs = reads
    
    for run in range(len(reads),1,-1):
        matches = {}
        for read1 in scs:
            for read2 in scs:
                if read1 != read2:
                    len_reads = min(len(read1),len(read2))
                    for i in range(len_reads,0,-1):
                        if read1[len(read1)-i:len(read1)] == read2[0:i]:
                            if i in matches:
                                matches[i][read1] = read2
                            else:
                                matches[i] = {read1:read2}
                            break
        for length in sorted(matches.keys(), reverse=True):
            for match1 in matches[length]:
                match2 = matches[length][match1]
                new = match1 + match2[length:len(match2)+1]
                scs.append(new)
                scs.remove(match1)
                scs.remove(match2)
                #print(match1,"+",match2,"-->",new)
                break
            break
    scs = scs[0]
    return scs

def main(filename):
    # read the sequences from the file
    sequences = read_sequences(filename)
    print("Read the sequences", file=stderr)

    # calculate the shortest common superstring
    superstring = calculate_scs(sequences)

    # print the result
    print(superstring)

if __name__ == "__main__":
    try:
        opts, args = getopt(argv[1:], "h", ["help"])
    except GetoptError as err:
        print(err)
        print(__doc__, file=stderr)
        exit(1) 

    for o, a in opts:
        if o in ("-h", "--help"):
            print(__doc__, file=stderr)
            exit()
        else:
            assert False, "unhandled option"

    if len(args) != 1:
        print(__doc__, file=stderr)
        exit(2)

    main(args[0])
