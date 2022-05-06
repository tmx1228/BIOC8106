#!/usr/bin/env python

"""
    usage:
        remove_transitive_edges [options] graph.dot
    where the options are:
        -h,--help : print usage and quit

    graph.dot is a file with a directed graph. The first row is always
        
        Digraph G {

    and the last row is always 

        }

    Each row in the middle section describes an edge. For example,

        1 -> 2

    is a directed edge from a node with the label '1' to the node with the label    '2'.
"""

from sys import argv, stderr
from getopt import getopt, GetoptError
from copy import deepcopy
from graph import *

def simplify(G):
    """Simplify the graph S by removing the transitively-inferrible edges.

    S is just a copy of G, which is the input to the graph. 
    """
    S = deepcopy(G)
    edges = {}
    nodes = G.nodes()
    for node1 in nodes:
        for node2 in nodes:
            if G.has_edge(node1,node2):
                if not node1 in edges:
                    edges[node1] = {node2:1}
                elif not node2 in edges[node1]:
                    edges[node1][node2] = 1
                else:
                    edges[node1][node2] = edges[node1][node2]+1
    similar_edges = deepcopy(edges)
    for node1 in edges:
        for node2 in edges[node1]:
            for node3 in edges[node1]:
                if G.has_edge(node2,node3):
                    if node3 in edges[node2]:
                        similar_edges[node1].pop(node3,None)
                    elif node2 in edges[node3]:
                        similar_edges[node1].pop(node2,None)   
    
    S = Graph(similar_edges)
    return S   

def main(filename):
    # read the graph from the input file
    graph = Graph(filename)
    print(f"Read the graph from {filename}", file=stderr)

    # simplify the graph by removing the transitively-inferrible edges
    simplified = simplify(graph)
    print(f"Simplified the graph", file=stderr)

    # print the simplified graph in the same format as the input file
    print(simplified)

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
