import networkx as nx
import random as r
import numpy as np
import matplotlib.pyplot as plt


t = 100     # length of the simulation --> must be at least 12 months for accuracy


def build_ba_model(n, m):
    # a function that builds a power-law graph (BA = scale-free)
    # n is the number of people, m is the number of connections of each node
    G = nx.generators.barabasi_albert_graph(n, m, seed=None)
    return G


# need to make sure that at least 1 person actually starts with TB !!
def add_attributes(G):
    # a function that adds initial attributes to each node
    for node in G:
        # 0 is susceptible, 1 is infected, 2 is recovered -> initialises a grid of susceptible nodes

        # start with all agents alive:
        G.nodes[node]['alive'] = [1]

        # chooses one of 0, 1, 2, with probs in p:
        #G.nodes[node]['state'] = list(np.random.choice([0, 1, 2], 1, p=[0.98734, 6.86e-3, 5.8e-3])) # check stats w PHRU
        G.nodes[node]['state'] = list(np.random.choice([0, 1, 2], 1, p=[0.9, 0.04, 0.06]))

        # assigns HIV+ status w prob p:
        G.nodes[node]['HIV_status'] = np.random.binomial(1, 0.013) # HIV incidence SA

        ## NEED TO WRITE FUNCTIONS TO ASSIGN THESE (TARGETED)
        # no one starts on TBT:
        G.nodes[node]['TB_treatment'] = [0]
        # no one starts on IPT:
        G.nodes[node]['IPT'] = [0]

        #need to assign TB preventative treatment (to 0 and 2) and TB treatment (only to 1)

    return(G)


def add_weights(G):
    # a function that adds weights to the edges of a graph
    for edge in G.edges:
        # assigns weights randomly with equal prob
        G[edge[0]][edge[1]]['weight'] = r.uniform(0, 1)
    return


def add_edge(lst, edge): 
    # function that adds all the visited edges to a list
    lst.append(edge)
    new_edge = [edge[1], edge[0]]
    lst.append(new_edge)
    return lst


def setup(n, m):
    G = build_ba_model(n, m)
    add_attributes(G)
    add_weights(G)
    #random_choose(G, 0.2) # percentages of edges that are going to be removed
    return G


## code for spreading of the infection

def infect(G, time, root, num, visited, weight):
    # a function that recursively visits all nodes and updates their states based on the states of their neighbours
        state(G, num, root, time, weight)
        for connection in G.edges(root):
            if connection not in visited:
                weight = G[connection[0]][connection[1]]['weight']
                visited = add_edge(visited, connection)  # adds the edge to a list of visited edges
                connection = connection[-1]  # gets the connected node
                infect(G, time, connection, G.nodes[root]['state'][time], visited, weight)
            else:
                return

# probabilities code:

# what needs to happen:
## if node is 0 = susceptible
## --> need to check if HIV+ or not
## --> if yes, prob of getting infected =?
## --> if no, prob of getting infected =?
# if node is 1 = TB disease
# --> need to check if on treatment or not. then in both cases:
# --> prob of death = ? (remove node from graph)
# --> prob of recovery (i.e. move to 2 (infected)) = ?
# if node is 2 = TB infection
# if on preventative treatment or not
# --> prob of re-infection (with disease)
# --> prob of death (?)
# IS IT POSSIBLE TO FULLY RECOVER FROM TB ??

def state(G, connected_state, node, time, weight): 
    # a function that determines the next state of the given node
    lst = G.nodes[node]['state'] # list of past states of the specific node for all times that have passed
    if lst[time] == 0:  # if the node is currently susceptible

        if len(lst)<(time+2): # if it hasn't already been visited

            # HIV_status
            if G.nodes[node]['HIV_status'] == 1: # HIV_pos
                if connected_state == 1: # if the connected node has TB disease
                  disease = np.random.binomial(1, 0.04*20*weight) # prob of infection with TB if HIV+
                  lst.append(disease)
                else:
                  lst.append(0)
            else: # HIV_neg
                if connected_state == 1:
                  disease = np.random.binomial(1, 0.04*weight) # prob of infection with TB if HIV-
                  lst.append(disease)
                else:
                  lst.append(0)

        else:
            if lst[time+1] == 0 and connected_state == 1: # if the node has been in contact with infected individuals
                # HIV_status
                if G.nodes[node]['HIV_status'] == 1:  # HIV_pos
                    if connected_state == 1:  # if the connected node has TB disease
                        disease = np.random.binomial(1, 0.04*20*weight)  # prob of infection with TB if HIV+
                        lst.append(disease)
                    else:
                        lst.append(0)
                else:  # HIV_neg
                    if connected_state == 1:
                        disease = np.random.binomial(1, 0.04*weight)  # prob of infection with TB if HIV-
                        lst.append(disease)
                    else:
                        lst.append(0)

    elif lst[time] == 1: # if agent has TB disease
        if len(lst) < (time+2):

            if G.nodes[node]['TB_treatment'] == 1:
                death = np.random.binomial(1, 0.03*weight)   # probability of death if TB disease and on treatment
                if death == 1:
                    weight_to_zero(G, node)
                    G.node[node]['alive'] = [0]
                elif recover(lst, node, G):
                    lst.append(2)
                else:
                    lst.append(1)
            else: # not on TBT
                death = np.random.binomial(1, 0.08*weight)  # probability of death if TB disease and not on treatment
                if death == 1:
                    weight_to_zero(G, node)
                elif recover(lst, node, G): # need to change recover function!!!
                    lst.append(2)
                else:
                    lst.append(1)

    elif lst[time] == 2:  # if agent is infected with TB
         if len(lst) < (time + 2):

             if G.nodes[node]['IPT'] == 1:
                 reinfect = np.random.binomial(1, 0.9*weight) + 1  # (1 - prob of re-infection with TB disease)
                 lst.append(reinfect)
             else:
                 lst.append(2)
    return


def spread(G, n):
    # function to spread a disease over a given number of time steps
    time = 0
    root = r.randint(0, n-1)     # chooses a random root at which to start checking the spread
    weight = 1  # initial weight for root infected node, not important
    while time < t:
        ## HERE IS WHERE WE WOULD LATER IMPLEMENT INTERVENTIONS:
        #if time == 0:
            #rm_edge_weight(G, 0.5) # remove edges keeping ones with highest weights only
            #random_choose(G, 0.2) # remove edges randomly - control
            #percolate(G) # remove edges based on edge betweeness centrality
        visited =[]
        state(G, G.nodes[root]['state'][time], root, time, weight)  # updates root
        infect(G, time, root, G.nodes[root]['state'][time], visited, weight)
        time += 1


# change this function to make sense for TB recovery
def recover(lst, node, G): 
    # function that determines whether or not an infected node recovers

    if G.nodes[node]['TB_treatment'] == 1:
        num = 2 # infectious period approx 2 months
    else:
        num = 24 # infectious period 24 months if no treatment

    summ = 0
    for i in lst:
        summ += i
    if summ > num:
        return True
    else:
        return False


# Data code 

def data(G):
    DATA = nx.get_node_attributes(G,'state') # dictionary -> node : [list of states of that node for all times]
    tot_S = [] # 0 = susceptible
    tot_TBd = [] # 1 = TB disease
    tot_TBi = [] # 2 = TB infection
    #tot_deaths = []
    data = [tot_S, tot_TBd, tot_TBi]
    
    for time in range(t):
        
        my_list = [elem[time] for elem in DATA.values()]    # list of states at a given time step
        
        s = my_list.count(0) # counts how many values are in state 0 - how many nodes are infected at time
        tot_S.append(s)
        
        i = my_list.count(1)
        tot_TBd.append(i)
        
        r = my_list.count(2)
        tot_TBi.append(r)
        
    print('peak', max(tot_TBd)) # print all necessary information
    print('occurs at', tot_TBd.index(max(tot_TBd)))
    print('never infected' , min(tot_S))
    print('no new infection' , tot_S.index(min(tot_S)))
    print('total that get it', max(tot_TBi))
    print('all recover at', tot_TBi.index(max(tot_TBi)))
    return data    

# percolation code
def weight_to_zero(G, node):
    # takes in a node and makes all its edge weights = 0: 
    for edge in G.edges(node): 
        G[edge[0]][edge[1]]['weight'] == 0 
    

def rm_edge_weight(G, threshold):
    # negates edges below a certain threshold (discard if their weight is below the threshold)
    for node in G:
        for edge in G.edges(node):
            if G[edge[0]][edge[1]]['weight'] <= threshold:
                G[edge[0]][edge[1]]['weight'] = 0
        
def percolate(G):
    # the function will remove edges from the network how to do it at some time t??
    remove = [] # edges to be removed
    centralities = nx.edge_betweenness_centrality(G)
    for edge in G.edges(): 
        if centralities[edge] > 0.1:
            remove.append(edge)
    for edge in remove: 
        G[edge[0]][edge[1]]['weight'] = 0 
    return 

def random_choose(G, cut_prob):
    # randomly negates certain edges of a graph 
    chosen = []
    for node in G: 
        for edge in G.edges(node): 
            x = np.random.binomial(1, cut_prob) # percentage of edges that will be cut in the beginning
            if x == 1:
                chosen.append(edge)
                G[edge[0]][edge[1]]['weight'] = 0 # USE THIS - SETTING EDGE WEIGHT TO ZERO TO AVOID DATA FUNTION NON WORKING
    #G.remove_edges_from(chosen)
    #print(chosen)
    return


def colour(G):
    colour_map = []
    for node in G:
        if G.nodes[node]['alive'] == 0:
          colour_map.append('gray')
        elif G.nodes[node]['state'] == 0:
          colour_map.append('green')
        elif G.nodes[node]['state'] == 1:
          colour_map.append('red')
        elif G.nodes[node]['state'] == 2:
          colour_map.append('orange')
    return colour_map


if __name__ == "__main__":
    n = 100  # S0 -> number of nodes at time t = 0.
    m = 3
    G = setup(n, m)

    colour(G)
    print(colour(G))

    spread(G,n)
    nx.draw(G) #node_color=colour_map) # --> draws the graph G
    plt.show()
    plt.plot(data(G)[0]) #plotting time evolution of nuber of susceptible
    plt.plot(data(G)[1]) #plotting time evolution of number of infected
    plt.plot(data(G)[2]) #plotting time evoluition of number of recovered
    plt.show()





# Questions for Workshops
# graph using a power-law distrib
# weights drawn from a uniform distribution — would this work? or is there a better way to find weights
# (can't delete nodes)
# how to assign the interventions?
# how to add deaths?


# For Research:
# stats for initialisation:
# - need % of pop initially infected (1), susceptible (0), recovered (2)
# - % of HIV+ people (could add a scale instead of just +/-)


# interventions:
# target most-connected nodes
# high-strength connections

# add prob of dying? for everyone?

# need to make sure that at least 1 person actually starts with TB !!