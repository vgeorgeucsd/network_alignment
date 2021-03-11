import networkx as nx
import matplotlib.pyplot as plt 
import re
import os

def make_prime(fname):
    G = nx.DiGraph()
    G=nx.read_edgelist(fname,create_using=nx.DiGraph())
    H = nx.Graph()
    for e in G.edges():
        src=e[0]
        tar=e[1]
        prime=e[1]+'p'
        H.add_edge(src, prime)
        H.add_edge(prime,tar)
    out_fname = os.path.splitext(fname)[0] + '_p' + '.el'
    print(out_fname)
    nx.write_edgelist(H, out_fname)
    
fname_directed = '/home/vikash/ssd2/SecondDrive/randomized_edge_lists/well0_day100000.txt'
make_prime(fname_directed)