import os
from math import log as log
import pandas as pd
import csv
import time
import numpy as np
import networkx as nx

cwd = os.getcwd()

def make_Orbit_Degree_Matrix_Dictionary(inpath):
    
    ''' Function Purpose: Creates a dictionary called the 'orbit_degree_matrices_dictionary'  whose keys are 
    'days' and values are the pair of file paths for the 'signature' file and node 'dictionary' file. 
    
    Functions Called:
    A. appendPath'''
    
    orbit_degree_matrices_dictionary = {}
    days = []
    
    for txt in os.listdir(inpath):
        # populate the graphlet_vector_matrix_dictionary with paths for the signatures and dictionaries
        
        ii = None
        
        if txt.startswith("well%s" %(well_number)):
            if txt.endswith(".signatures.txt"):
                ii = 0
            elif txt.endswith(".dictionary.txt"):
                ii = 1
            else:
                continue
            
            orbit_degree_matrices_dictionary, days = appendPath(orbit_degree_matrices_dictionary, days, txt, inpath, ii)
    
    #graphlet_vector_matrix_dictionary = dict(sorted(graphlet_vector_matrix_dictionary.items(), key=lambda x: x[0]))

    return orbit_degree_matrices_dictionary, days

def appendPath(orbit_degree_matrices_dictionary, days, txt, inpath, ii):
    ''' This appends the file path to either the signatures set (ie graphlet_degree_vectors) 
    or the the dictionary set'''
    
    day = int(txt.split('.')[0].split('_')[1].split('y')[1])
    
    if not day in days:
        
        days.append(day)
        
        days = sorted(days) # 'days' stores integer values in increasing sequence
        
        orbit_degree_matrices_dictionary[day] = [None, None]
    
    path = inpath + txt
    
    orbit_degree_matrices_dictionary[day][ii] = path
    
    return orbit_degree_matrices_dictionary, days
    
def store_Data_For_TimeGraphs(orbit_degree_matrices_dictionary, current_day, next_day, d, gap, orbit_degree_matrix_1, orbit_degree_matrix_2, node_labels_1, node_labels_2):
    
    ''' Function Purpose: To extract the signature files of two graphs with time points.
    
        Variables: 
        1. graphlet_vector_matrix_dictionary := a dictionary which holds information about each node's orbit degree vector with the corresponding node label stored in the dictionary. 
        2. days := this is an ordered array with all the time stamps (unit: day) 
        3. d    := the 'days' array index. 
        4. gap  := this value tells us how far apart (units: day) to sample the graphs. Ex: 0 means adjacent time stamps.
        5. orbit_degree_matrix_1   := for graph #1, this is its orbit degree matrix (whose rows correspond to a specific node and columns are the orbits). 
        6. orbit_degree_matrix_2   := for graph #2, this is its orbit degree matrix (...)
        7. node_labels_1 := for graph #1, this dictionary stores the node label for each of its rows. (Ie. The first row does not necessarily mean node 1)
        8. node_labels_2 := for graph #2, this dictionary stores its node labels. 
        
        Functions Called:
        A. read_OrbitDegree_Matrix
        B. read_Node_Labels           '''

    current_graph = orbit_degree_matrices_dictionary[current_day]
    next_graph = orbit_degree_matrices_dictionary[next_day]


    if gap == 0 and d > 0:
        # switch orbit_degree_matrix_2 to orbit_degree_matrix_1 ie use the info from graph i and compute the new info for graph i+1
        orbit_degree_matrix_1 = orbit_degree_matrix_2
        node_labels_1 = node_labels_2

        orbit_degree_matrix_2 = np.asarray(read_OrbitDegree_Matrix(next_graph[0]))
        node_labels_2 = np.asarray(read_Node_Labels(next_graph[1]))


    else:

        orbit_degree_matrix_1 = np.asarray(read_OrbitDegree_Matrix(current_graph[0]))
        node_labels_1 = np.asarray(read_Node_Labels(current_graph[1]))  #read_OrbitDegree_Matrices := it reads the orbit degree vectors of all the nodes in a graph

        orbit_degree_matrix_2 = np.asarray(read_OrbitDegree_Matrix(next_graph[0]))
        node_labels_2 = np.asarray(read_Node_Labels(next_graph[1]))


    out = outpath + 'day' + str(current_day) + '_' + 'day' + str(next_day) + '.csv'
    # IF MATRIX, CHANGE HERE. 

    return orbit_degree_matrix_1, orbit_degree_matrix_2, node_labels_1, node_labels_2, out

def read_OrbitDegree_Matrix(file):
    
    ''' Function Purpose: Given an orbit degree matrix file, this function reads it and stores it into a list.'''
    
    orbit_degree_matrix = []

    fRead = open(file, 'r')

    for line in fRead:
        splitted = line.strip().split(' ')
        orbit_degree_matrix.append([int(value) for value in splitted])
    fRead.close()

    return orbit_degree_matrix

def read_Node_Labels(file):
    
    ''' Function Purpose: This functions outputs an array of all the node labels for the rows in the orbit degree matrix'''
    
    node_labels = []

    fRead = open(file, 'r')

    for line in fRead:
        line = line.splitlines()[0].split(' ')
        index = int(line[0])
        node = int(line[1])
        node_labels.append(node) # We subtract one to counter the original offset of +1. Another way to do this is to rerun everything with 'ucsd_txt' instead of 'ucsd_txt_nonzero'
    fRead.close()
    
    return node_labels

def populateSimilarityMatrix(orbit_degree_matrix_1, orbit_degree_matrix_2, SimMatrix):

    ''' Function Purpose: Calculates the similarity between two orbit degree matrices. Stores similarity values in an empty matrix called SimMatrix 
    
        Variables: 
        1. orbit_degree_matrix_1 := the orbit degree matrix for graph 1
        2. orbit_degree_matrix_2 := the orbit degree matrix for graph 2
        3. SimMatrix := an empty matrix for the similarity values.       
        
        Functions Called:
        A. Average_OrbitDistance'''
    
    for node1 in range(0,len(orbit_degree_matrix_1)): #loops over all orbit degree vectors in graph 1
        
        for node2 in range(0,len(orbit_degree_matrix_2)): # loops over all orbit degree vectors in graph 2
            
            
            DGDVS_value = Average_orbitDistance(orbit_degree_matrix_1, orbit_degree_matrix_2, node1, node2)
                # compute the Directed Graphlet Degree Vector Similarity
            
            SimMatrix[node1, node2] = DGDVS_value 
                # store the similarity value in the matrix entry for node 1 and node 2    
                
            ###DiffMatrix[node1, node2] = DGDVS_value
                ### store the DISsimilarity value in the matrix entry for node 1 and node 2
                
    return SimMatrix

def Average_orbitDistance(orbit_degree_matrix_1, orbit_degree_matrix_2, node1, node2):
    
    ''' Function Purpose: Computes similarity between to orbit degree vectors based on the average over all orbit distances.
        
        Variables:
        1. orbit_degree_matrix_1 := orbit degree matrix for graph 1
        2. orbit_degree_matrix_2 := orbit degree matrix for graph 2
        3. node1 := Node from orbit degree matrix 1
        4. node2 := Node from orbit degree matrix 2  
        
        Functions Called:
        A. orbitDistance       '''
    
    SumRow = 0
    SumWi = 0
    
    for orbit in range(0, 128):
    # 128 orbits to be considered.
    
        orbit_distance, wi = orbitDistance(orbit_degree_matrix_1, orbit_degree_matrix_2, node1, node2, orbit)
        SumRow = SumRow + orbit_distance
        SumWi = SumWi + wi
    
    Distance = SumRow/SumWi
    GDS_value = 1 - Distance
    # This is the distance between node 1 and node 2 with respect to their graphlet degree vectors
    
    return GDS_value

def orbitDistance(orbit_degree_matrix_1, orbit_degree_matrix_2, node1, node2, orbit):
    
    ''' Function Purpose: Computes the difference between two nodes with respect to a single orbit.
                          Takes into account the absolute difference between orbit degrees. 
    
        Variables: 
        1. orbit_degree_matrix_1 := orbit degree matrix for graph 1
        2. orbit_degree_matrix_2 := orbit degree matrix for graph 2
        3. node1 := Node from orbit degree matrix 1
        4. node2 := Node from orbit degree matrix 2 
        5. orbit := The specific orbit to be used.   
        
        Functions Called: 
        A. calcWeight
        B. calcOi        '''
    
    
    wi = calcWeight(calcOi(orbit))
    
    num = abs(orbit_degree_matrix_1[node1, orbit] - orbit_degree_matrix_2[node2, orbit])
    denom = max(orbit_degree_matrix_1[node1, orbit], orbit_degree_matrix_2[node2, orbit])
    
    if denom == 0:
        orbit_distance = 0
    else:
        orbit_distance = wi*num/denom
    
    return orbit_distance, wi

def calcWeight(oi):
    wi = 1 - (log(oi))/(log(129))
    return wi

def calcOi(orbit):
    
    ''' Function Purpose: Receives an orbit number and outputs the number of other orbits in which it depends on. 
        Ex: A center node (/orbit) in a star depends on the number of times the node is a source node in a directed edge.'''
    oi_1 = [0, 1]
    oi_2 = [2,4,5,6,7,8,10,11]
    oi_3 = [3,9,12,13,16,17,20,21,24,25,28,29,30,31,32,33,36,60,64,68,72,118,120,125,128]
    oi_4 = [19,22,23,27,35,38,44,45,47,48,52,56,57,61,63,66,67,70,73,74,76,77,78,80,85,87,91,92,95,96,100,105,113,115]
    oi_5 = [14,15,18,26,34,37,42,43,49,50,53,54,58,59,62,65,69,71,88,90,93,94,99,103,107,111,121,122,126,127]
    oi_6 = [39,40,41,46,51,55,75,79,81,82,83,84,86,89,97,98,101,102,104,106,108,109,110,112,114,116,117,119,123,124]
    if orbit in oi_1:
        return 1
    elif orbit in oi_2:
        return 2
    elif orbit in oi_3:
        return 3
    elif orbit in oi_4:
        return 4
    elif orbit in oi_5:
        return 5
    elif orbit in oi_6:
        return 6
    else:
        return 0
    
def convertSimMatrix2SimList(matrix):
    
    ''' Function Purpose: Converts the similarity matrix pandas dataframe and converts it to a list of similarities
    where the first column is the node 1 from graph1, second column is the node 2 from graph2, and the third
    column is the similarity value. The order goes from all the pairings with node1 to all the other nodes in graph 2.'''
    
    sim_matrix = matrix
    
    node_pair_similarities = []
    
    graph1_node_array = list(sim_matrix.index)
    
    graph2_node_array = list(sim_matrix.columns)

    for i in range(0, len(graph1_node_array)): # index through the vertices of graph 1, (rows)

        node1 = graph1_node_array[i]

        for j in range(0, len(graph2_node_array)): # index through the vertices of graph 2, (columns)

            node2 = graph2_node_array[j]

            sim_value = sim_matrix.loc[node1, node2] # extract the entry value in (row, column)

            node_pair_similarities.append([node1, node2, sim_value]) # append to the list of nodes and their similarities
            
    return node_pair_similarities

def In_Out_Paths(GAP, OUTPUT_FORMAT, cwd):
    
    ''' Function Purpose: Picks the correct inpath and outpaths. The out paths considers if user wants gaps in days or not. 
    Ie to consider consecutive time graphs or some other extension of time. Also, considers if user wants lists or not.
    
    Varaibles:
    1. GAP := How far apart are the time graphs' time points?
    2. OUTPUT_FORMAT := Do you want it in a list format?
    3. cwd := current working directory '''
    
    inpath = None
    outpath = None
    outpath2 = None
    
    ## SET INPATHS & OUTPATHS 
    # Outpaths differentiates y/n gaps, and y/n lists
        
    ## Inpath:
    inpath = 'CHANGE_INPUT_PATH'

    ## Outpath:
    if GAP == 0:
        GDS_dir = '/GDS_days/'

    ## Outpath_gaps:
    elif GAP > 0:
        parent_dir = '/GDS_gap_days/'

        gap_dir = 'GDS_gap_' + str(GAP) + '_days/'
        
        GDS_dir = parent_dir + gap_dir                    # Specifies the gap in the directory
        
    ## Outpath needs gaps or not:
    outpath = 'CHANGE_OUTPUT_PATH' + GDS_dir + 'well' + well_number + '/'  # where to place the similarity matrices

    ## Asks if directory exists:
    if os.path.isdir(outpath) == False:
        os.makedirs(outpath)
        print("Directory '% s' created" % outpath)

    ## Outpath_lists:
    if OUTPUT_FORMAT == 'list':

        outpath2 = 'CHANGE_OUTPUT_PATH' + GDS_dir + 'well' + well_number + '/'  # where to place the similarity matrices

        if os.path.isdir(outpath2) == False:
            os.makedirs(outpath2)
            print("Directory '% s' created" % outpath2)
    
    
    return inpath, outpath, outpath2

WELLS = ['22', '23', '34']

## SET OPTIONS ------------
#for GAP in range(0, 4):
    
GAP = CHANGE_GAP  # CHANGE # Nonnegative number indicating the difference between time of compared graphs. '0' indicates adjacent time graphs.

OUTPUT_FORMAT = 'matrix'   # If 'list', then does the procedure to convert matrix to list. Right now, 'list' produces both matrix and list format.

STOP_WATCH = False       # If 'True', then does timer operations



#for w in WELLS:

well_number = 'CHANGE_WELL_NUMBER' # CHANGE


## SET INPATHS & OUTPATHS 

# Paths:
inpath, outpath, outpath2 = In_Out_Paths(GAP, OUTPUT_FORMAT, cwd)



## INITIALIZE PARAMETERS, ARRAYS, & OPTIONS

## Parameters:
gap = GAP

## Arrays:
orbit_degree_matrix_1 = []
orbit_degree_matrix_2 = []

node_labels_1 = []
node_labels_2 = []

## Options:
if OUTPUT_FORMAT == 'list':
    Similarity_Matrices_Set = []

if STOP_WATCH == True:
    array_time_for_matrix_df = []
    array_time_for_list_df = []
    array_time_for_matrix2list_conversion = []



## COMPUTING THE SIMILARITY BETWEEN TIME GRAPHS ACCORDING TO THEIR ORBIT DEGREE MATRICES

## Get all the orbit degree matrices for the well (from txt to dictionary):
orbit_degree_matrices_dictionary, days = make_Orbit_Degree_Matrix_Dictionary(inpath)


for d in range(0,len(days)-1): 

    if (d + gap + 1) < len(days):     # This ensures we do not index out of range.

        # Set days.
        current_day = days[d]
        next_day = days[d + gap + 1] # 'gap' = 0 refers to adjacent time points.

        ## Start the timer:
        if STOP_WATCH == True:
            start = time.time()

        ## Get the orbit degree matrices (from dictionary to arrays):
        orbit_degree_matrix_1, orbit_degree_matrix_2, node_labels_1, node_labels_2, out = store_Data_For_TimeGraphs(orbit_degree_matrices_dictionary, current_day, next_day, d, gap, orbit_degree_matrix_1, orbit_degree_matrix_2, node_labels_1, node_labels_2)


        #M = max(len(SigFile1), len(SigFile2)) # M := maximum node set size between the two graphs.      

        # DiffMatrix = np.ones((M , M))  
        # DiffMatrix has dimensions equal to the maximum node set. Fill entries with '1' since '1' indicates maximum distance/difference. This is mostly related to unaligned nodes. 
        # SimMatrix = np.zeros((M, M))

        ## Compute the matrix with similarity values:
        SimMatrix = np.zeros((len(orbit_degree_matrix_1), len(orbit_degree_matrix_2)))

        SimMatrix = populateSimilarityMatrix(orbit_degree_matrix_1, orbit_degree_matrix_2, SimMatrix)

        df = pd.DataFrame(data=SimMatrix, index=node_labels_1, columns=node_labels_2)

        df.sort_index(inplace=True)
        
        df.sort_index(axis=1,inplace=True)
        
        df.to_csv(out,header=False,index=False)

        ## Get time to build matrix:
        if STOP_WATCH == True:

            stop1 = time.time()



        ## Get similarity lists:
        if OUTPUT_FORMAT == 'list':

            #Similarity_Matrices_Set.append(df)

            out2 = outpath2 + 'day' + str(current_day) + '_' + 'day' + str(next_day) + '.csv' # CHANGE 

            node_pair_similarities = convertSimMatrix2SimList(df) 

            df_2 = pd.DataFrame(data=node_pair_similarities, index=None, columns=['node 1', 'node 2', 'similarity'])

            df_2.to_csv(out2)


        ## Get times:
        if STOP_WATCH == True:

            stop2 = time.time()

            time_for_matrix_df = stop1 - start

            if OUTPUT_FORMAT == 'list':

                time_for_list_df = stop2 - start

                time_for_matrix2list_conversion = stop2 - stop1

            #array_time_for_matrix_df.append(time_for_matrix_df)
            #array_time_for_list_df.append(time_for_list_df)
            #array_time_for_matrix2list_conversion.append(time_for_matrix2list_conversion)

            print('time_for_matrix_df: ',time_for_matrix_df, '\n', 'time_for_list_df: ', time_for_list_df, '\n', 'time_for_matrix2list_conversion', time_for_matrix2list_conversion)

node_labels_dictionary = None
    # (str) File path to get the dictionary. 

class OrbitDegreeSimilarityGraph:
    ''' 
    This takes the orbit degree matrices of two graphs and 
    compares them via their nodes' orbit degree vectors similarities.

    '''
    number_of_orbits = 129
    # (int) Total number of orbits
    
    
    
    def __init__(self, orbit_degree_matrix_1, orbit_degree_matrix_2, node_labels_dictionary):  
        
        self.orbit_degree_matrix_1 = orbit_degree_matrix_1
        self.orbit_degree_matrix_2 = orbit_degree_matrix_2
        self.orbits = range(0, number_of_orbits)
        self.similarity_graph = nx.complete_bipartite_graph(len(self.orbit_degree_matrix_1), (len(self.orbit_degree_matrix_2)))
        #self.node_labels_for_graph_1 = 
        # label the nodes.
        # get mapping. (dictionary to convert node labels)
        
        self.similarity_graph = nx.relabel_nodes(self.similarity_graph, node_labels, copy=False)
        
        
        
        self.similarity_list = np.zeros((len(self.orbit_degree_matrix_1) * len(self.orbit_degree_matrix_2)), 2)
        self.SimMatrix = np.zeros((len(self.orbit_degree_matrix_1), len(self.orbit_degree_matrix_2)))
        
    
    def makeSimilarityGraph(self):
        
        for v1 in range(0,len(self.orbit_degree_matrix_1)): 
        #loops over all orbit degree vectors in graph 1

            for v2 in range(0,len(self.orbit_degree_matrix_2)): 
            # loops over all orbit degree vectors in graph 2
                
                DGDVS_value = Average_orbitDistance(self.orbit_degree_matrix_1, self.orbit_degree_matrix_2, v1, v2)
                # compute the Directed Graphlet Degree Vector Similarity
                
                self.similarity_graph[v1][v2]['weight'] = DGDVS_value

    
    
    def populateSimilarityMatrix(self): 

        ''' Function Purpose: Calculates the similarity between two orbit degree matrices. Stores similarity values in an empty matrix called SimMatrix 
    
        Variables: 
        1. self.orbit_degree_matrix_1 := the orbit degree matrix for graph 1
        2. self.orbit_degree_matrix_2 := the orbit degree matrix for graph 2
        3. self.SimMatrix := an empty matrix for the similarity values.       
        
        Functions Called:
        A. Average_OrbitDistance'''

        for node1 in range(0,len(self.orbit_degree_matrix_1)): 
        #loops over all orbit degree vectors in graph 1

            for node2 in range(0,len(self.orbit_degree_matrix_2)): 
            # loops over all orbit degree vectors in graph 2

                DGDVS_value = Average_orbitDistance(self.orbit_degree_matrix_1, self.orbit_degree_matrix_2, node1, node2)
                # compute the Directed Graphlet Degree Vector Similarity

                self.SimMatrix[node1, node2] = DGDVS_value 
                # store the similarity value in the matrix entry for node 1 and node 2    

        
        
    def Average_orbitDistance(self, node1, node2):
    
        ''' Function Purpose: Computes similarity between to orbit degree vectors based on the average over all orbit distances.
        
        Variables:
        1. orbit_degree_matrix_1 := orbit degree matrix for graph 1
        2. orbit_degree_matrix_2 := orbit degree matrix for graph 2
        3. node1 := Node from orbit degree matrix 1
        4. node2 := Node from orbit degree matrix 2  
        
        Functions Called:
        A. orbitDistance       '''

        SumRow = 0
        SumWi = 0

        for orbit in self.orbits:

            orbit_distance, wi = self.orbitDistance(self.orbit_degree_matrix_1, self.orbit_degree_matrix_2, node1, node2, orbit)
            SumRow = SumRow + orbit_distance
            SumWi = SumWi + wi

        Distance = SumRow/SumWi
        GDS_value = 1 - Distance
        # This is the distance between node 1 and node 2 with respect to their graphlet degree vectors

        return GDS_value
        
        
    def orbitDistance(self, node1, node2, orbit):
    
        ''' Function Purpose: Computes the difference between two nodes with respect to a single orbit
    
        Variables: 
        1. orbit_degree_matrix_1 := orbit degree matrix for graph 1
        2. orbit_degree_matrix_2 := orbit degree matrix for graph 2
        3. node1 := Node from orbit degree matrix 1
        4. node2 := Node from orbit degree matrix 2 
        5. orbit := The specific orbit to be used.   
        
        Functions Called: 
        A. calcWeight
        B. calcOi        '''


        wi = self.calcWeight(self.calcOi(orbit))

        num = abs(log(orbit_degree_matrix_1[node1, orbit] + 1) - log(orbit_degree_matrix_2[node2, orbit] + 1))
        denom = log(max(orbit_degree_matrix_1[node1, orbit], orbit_degree_matrix_2[node2, orbit]) + 2)

        orbit_distance = wi*num/denom

        return orbit_distance, wi
 
    
    def calcOi(orbit):
    
        ''' Function Purpose: Receives an orbit number and outputs the number of other orbits in which it depends on. 
        Ex: A center node (orbit) in a star depends on the number of times the node is a source node in a directed edge.'''
        
        oi_1 = [0, 1]
        oi_2 = [2,4,5,6,7,8,10,11]
        oi_3 = [3,9,12,13,16,17,20,21,24,25,28,29,30,31,32,33,36,60,64,68,72,118,120,125,128]
        oi_4 = [19,22,23,27,35,38,44,45,47,48,52,56,57,61,63,66,67,70,73,74,76,77,78,80,85,87,91,92,95,96,100,105,113,115]
        oi_5 = [14,15,18,26,34,37,42,43,49,50,53,54,58,59,62,65,69,71,88,90,93,94,99,103,107,111,121,122,126,127]
        oi_6 = [39,40,41,46,51,55,75,79,81,82,83,84,86,89,97,98,101,102,104,106,108,109,110,112,114,116,117,119,123,124]
        
        if orbit in oi_1:
            return 1
        elif orbit in oi_2:
            return 2
        elif orbit in oi_3:
            return 3
        elif orbit in oi_4:
            return 4
        elif orbit in oi_5:
            return 5
        elif orbit in oi_6:
            return 6
        else:
            return 0
        
       
    def calcWeight(orbit):
        
        ''' Function Purpose: Computes the orbit weight.'''
        
        wi = 1 - (log(orbit))/(log(number_of_orbits))
        return wi