import networkx as nx
import random as r
import numpy as np
import matplotlib.pyplot as plt

t = 12*5     # length of the simulation --> must be at least 12 months for accuracy of stats and probs

## initialising statistics (percentages):
p_latentTB = 0.148              # unpub source
p_TBD = 0.041                   # unpub source
#p_TBD = 3.13e-3                # https://www.atsjournals.org/doi/full/10.1164/rccm.200409-1200OC
p_susceptible = 1 - p_latentTB - p_TBD
p_hiv_pos = 0.144               # unpub source

## probabilities:
#p_TBD_if_hiv+ = 0.106          # b
p_TBD_close_contacts = 0.04     # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6534268/
p_TBD_mortality = 0.06          # http://www.samj.org.za/index.php/samj/article/view/12734
p_TBD_death_on_treatment = 0.05 # conf by Dr
p_TBD_death_no_treatment = 0.2  # conf by Dr (hosp deaths) (20%)
p_reinfection_TBD = 0.18        # https://www.atsjournals.org/doi/full/10.1164/rccm.200409-1200OC

## treatments:
IPT_efficacy = 0.55             # https://journals.lww.com/aidsonline/fulltext/2003/09260/efficacy_of_secondary_isoniazid_preventive_therapy.7.aspx
death_hiv_neg_no_treatment = 0.45 # WHO: https://www.who.int/news-room/fact-sheets/detail/tuberculosis
death_hiv_pos_no_treatment = 0.9  # WHO: https://www.who.int/news-room/fact-sheets/detail/tuberculosis
p_self_cure_hiv_neg_no_TBT = 1/3  # conf by Dr


def build_ba_model(n, m):
    # a function that builds a power-law graph (BA = scale-free)
    # n is the number of people, m is the average number of connections of each node
    G = nx.generators.barabasi_albert_graph(n, m, seed=None)
    return G


# need to make sure that at least 1 person actually starts with TB!
def add_attributes(G):
    # a function that adds initial attributes to each node
    # lists so time evolution can be added
    for node in G:

        # start with all agents alive:
        G.nodes[node]['alive'] = [1] # make this a list over the times

        # 0 is susceptible, 1 is TBD, 2 is TB infected
        # chooses one of 0, 1, 2, with probs in p:
        #G.nodes[node]['state'] = list(np.random.choice([0, 1, 2], 1, p=[0.98734, 6.86e-3, 5.8e-3])) # check stats w PHRU
        G.nodes[node]['state'] = list(np.random.choice([0, 1, 2], 1, p=[p_susceptible, p_TBD, p_latentTB]))

        # assigns HIV+ status w prob p:
        G.nodes[node]['HIV_status'] = [np.random.binomial(1, p_hiv_pos)] # HIV incidence SA

        ## Interventions
        # no one starts on TBT:
        G.nodes[node]['TB_treatment'] = [0]
        # no one starts on IPT:
        G.nodes[node]['IPT'] = [0]
        # need to assign TB preventative treatment (to 0 and 2) and TB treatment (only to 1)

    return(G)


def add_weights(G):
    # a function that adds weights to the edges of a graph
    centralities = nx.edge_betweenness_centrality(G) # list of btwness centralities of edges
    for edge in G.edges:
        # assigns weights according to 'preferential weighting' (edges in clusters have greater weighting)
        G[edge[0]][edge[1]]['weight'] = centralities[edge]*10
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
    return G


## code for spreading the infection

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

def state(G, connected_state, node, time, weight): 
    # a function that determines the next state of the given node
    lst = G.nodes[node]['state'] # list of past states of the specific node for all times that have passed
    if lst[time] == 0:  # if the node is currently susceptible

        if len(lst)<(time+2): # if it hasn't already been visited

            # HIV_status
            if G.nodes[node]['HIV_status'] == [1]: # HIV_pos
                if connected_state == 1: # if the connected node has TB disease
                  disease = np.random.binomial(1, p_TBD_close_contacts*20*weight) # prob of infection with TB if HIV+
                  lst.append(disease)
                else:
                  lst.append(0)
            else: # HIV_neg
                if connected_state == 1:
                  disease = np.random.binomial(1, p_TBD_close_contacts*weight) # prob of infection with TB if HIV-
                  lst.append(disease)
                else:
                  lst.append(0)

        else:
            if lst[time+1] == 0 and connected_state == 1: # if the node has been in contact with infected individuals
                # HIV_status
                if G.nodes[node]['HIV_status'] == [1]:  # HIV_pos
                    if connected_state == 1:  # if the connected node has TB disease
                        disease = np.random.binomial(1, p_TBD_close_contacts*20*weight)  # prob of infection with TB if HIV+
                        lst.append(disease)
                    else:
                        lst.append(0)
                else:  # HIV_neg
                    if connected_state == 1:
                        disease = np.random.binomial(1, p_TBD_close_contacts*weight)  # prob of infection with TB if HIV-
                        lst.append(disease)
                    else:
                        lst.append(0)

    elif lst[time] == 1: # if agent has TB disease
        if len(lst) < (time+2):

            if G.nodes[node]['TB_treatment'] == [1]:
                death = np.random.binomial(1, p_TBD_death_on_treatment)
                if death == 1 or type(recover(lst, node, G)) == str: #includes risks of death based on HIV status
                    weight_to_zero(G, node)
                    G.nodes[node]['alive'] = [0]
                    lst.append(1)
                elif recover(lst, node, G):
                    lst.append(2)
                    # give treatment
                else:
                    lst.append(1)
            else: # not on TBT
                death = np.random.binomial(1, p_TBD_death_no_treatment)
                if death == 1 or type(recover(lst, node, G)) == str:
                    weight_to_zero(G, node)
                    G.nodes[node]['alive'] = [0]
                    lst.append(1)
                elif recover(lst, node, G):
                    lst.append(2)
                else:
                    lst.append(1)

    elif lst[time] == 2:  # if agent is infected with TB
         if len(lst) < (time + 2):

             if G.nodes[node]['IPT'] == [1]: # on IPT
                 if connected_state == 1:  # if the connected node has TB disease
                     reinfect = np.random.binomial(1, p_reinfection_TBD*(1-IPT_efficacy)*weight) + 1
                     lst.append(reinfect)
                 else:
                     lst.append(2)
             else:  # not on IPT
                 if connected_state == 1:
                     reinfect = np.random.binomial(1, p_reinfection_TBD*weight) + 1
                     lst.append(reinfect)
                 else:
                     lst.append(2)

         else:
             if lst[time + 1] == 0 and connected_state == 1:  # if the node has been in contact with infected individuals
                 if G.nodes[node]['IPT'] == [1]:  # on IPT
                     if connected_state == 1:  # if the connected node has TB disease
                         reinfect = np.random.binomial(1, p_reinfection_TBD*(1-IPT_efficacy)*weight) + 1
                         lst.append(reinfect)
                     else:
                         lst.append(2)
                 else:  # not on IPT
                     if connected_state == 1:
                         reinfect = np.random.binomial(1, p_reinfection_TBD*weight) + 1  # ???
                         lst.append(reinfect)
                     else:
                         lst.append(2)

    return


def spread(G, n, lst):
    # function to spread a disease over a given number of time steps
    time = 0
    root = r.randint(0, n-1)     # chooses a random root at which to start checking the spread
    weight = 1  # initial weight for root infected node, not important

    time_imp = lst[0]  # time
    type = lst[1]      # type
    letter = lst[2]    # letter
    while time < t:

        ## IMPLEMENT INTERVENTIONS at certain time steps:
        if time == time_imp:
            if type == 'wide':
                wide_treatment(G, letter, time)
                print('wide')
            elif type == 'cluster':
                cluster_treatment(G, letter, time)

        visited =[]
        state(G, G.nodes[root]['state'][time], root, time, weight)  # updates root
        infect(G, time, root, G.nodes[root]['state'][time], visited, weight)
        time += 1


## INTERVENTIONS

def wide_treatment(G, str, time):
    # function to give treatment to everyone who needs it
    # in a sense, implies limitless resources (nurses, time, follow-up, money etc)
    # str is a letter (in A, B, C) indicating the type of treatment:
        # str == A --> plan A: implement both IPT and TBT
        # str == B --> plan B: implement only TBT
        # str == C --> plan C: implement only IPT

    for node in G.nodes:
        # give all agents with TB disease treatment, make sure they aren't on IPT:
        if G.nodes[node]['state'][time] == 1 and (str == 'A' or str == 'B'):
            G.nodes[node]['TB_treatment'] = [1]
            G.nodes[node]['IPT'] = [0]
        # give all agents with TB infection OR HIV+ IPT if not on TBT
        elif (G.nodes[node]['state'][time] == 2 or G.nodes[node]['HIV_status'] == [1]) and G.nodes[node]['TB_treatment'] == [0] and (str == 'A' or str == 'C'):
            G.nodes[node]['IPT'] = [1]
    return

def cluster_treatment(G, str, time):
    # function to give treatment only to clusters
    # what does this mean realistically?
    # num as as before
    for node in G.nodes:
        if nx.clustering(G)[node] > 0.1 and G.nodes[node]['state'][time] == 1 and (str == 'A' or str == 'B'):
            G.nodes[node]['TB_treatment'] = [1]
            G.nodes[node]['IPT'] = [0]
        elif nx.clustering(G)[node] > 0.1 and (G.nodes[node]['state'][time] == 2 or G.nodes[node]['HIV_status'] == [1]) and G.nodes[node]['TB_treatment'] == [0] and (str == 'A' or str == 'C'):
            G.nodes[node]['IPT'] = [1]
    return


def vaccinate():
# decrease tie weights
# vaccination not included because its effects are negligible particularly in the adult population
    pass


def recover(lst, node, G): 
    # function that determines whether or not an infected node recovers
    summ = 0
    for i in lst:
        summ += i

    if G.nodes[node]['TB_treatment'] == [1]:
        num = 2  # infectious period approx 2 months (make this a longer choice?)
        if summ > num:
            return True # always recover after 2 month of treatment
        else:
            return False

    # if HIV infected, stay TBD or die (eventually) with no treatment
    elif G.nodes[node]['HIV_status'] == [1]:
        death = np.random.binomial(1, death_hiv_pos_no_treatment) # probability of dying from TB if HIV+ and NOT on treatment (check)
        if death == 1:
            return 'death'
        else:
            return False

    # HIV seronegative and TBD with no TB treatment - die/self cure/continue with HIV w prob 1/3 each
    else:
        num = 24 # infectious period 24 months if no treatment and HIV+
        if summ > num:
            death = np.random.binomial(1, death_hiv_neg_no_treatment) # or 1/3?
            self_cure = np.random.binomial(1, 1 / 3)  # from epi
            if death == 1:
                return 'death'
            elif self_cure == 1:
                return True
            else:
                return False


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

# Data code 

def data(G):
    DATA = nx.get_node_attributes(G,'state') # dictionary -> node : [list of states of that node for all times]
    tot_S = [] # 0 = susceptible
    tot_TBd = [] # 1 = TB disease
    tot_TBi = [] # 2 = TB infection

    DEATHS = nx.get_node_attributes(G,'alive')
    tot_deaths = []

    data = [tot_S, tot_TBd, tot_TBi, tot_deaths]
    
    for time in range(t):
        
        my_list = [elem[time] for elem in DATA.values()]    # list of states at a given time step
        
        s = my_list.count(0) # counts how many values are in state 0 - how many nodes are infected at time
        tot_S.append(s)
        
        i = my_list.count(1)
        tot_TBd.append(i)
        
        r = my_list.count(2)
        tot_TBi.append(r)

    death_list = [elem for elem in DEATHS.values()]
    d = death_list.count([0])
    tot_deaths.append(d)

    # changing from 0->1 or 2->1:
    move_to_TBD = 0
    move_to_LTB = 0
    TBD_recoveries = 0
    LTB_recoveries = 0
    for node in DATA:
        for i in range(len(DATA[node])):
            if i < len(DATA[node]) - 1:
                # if person gets active TBD
                if (DATA[node][i] == 0 and DATA[node][i+1] == 1) or (DATA[node][i] == 2 and DATA[node][i+1] == 1):
                    move_to_TBD += 1
                # if person becomes infected with LTB
                #elif DATA[node][i] == 0 and DATA[node][i+1] == 2:
                    #move_to_LTB += 1
                    #break # so that the same person moving is not counted more than once
                elif DATA[node][i] == 1 and DATA[node][i+1] == 2:
                    TBD_recoveries += 1
                #elif DATA[node][i] == 2 and DATA[node][i+1] == 0:
                    #LTB_recoveries += 1

    # print all necessary information
    print('# deaths attributable to TBD:', tot_deaths[-1])
    print('# people who contracted active TBD:', move_to_TBD)
    #print('# people who were infected with LTB:', move_to_LTB)
    print('# recoveries from TBD', TBD_recoveries)
    #print('# recoveries from LTB', LTB_recoveries)

    print('# infectious people at end of trial:', tot_TBd[-1])
    print('Maximum # of people with TBD at a point:', max(tot_TBd), 'occurs at', tot_TBd.index(max(tot_TBd)), 'months')
    print('Maximum # of people with TBI at a point:', max(tot_TBi), 'occurs at', tot_TBi.index(max(tot_TBi)), 'months')
    print('# of people never infected', tot_S[-1])

    print('-----------------------------------------------------')

    return data


def colour(G, t):
    colour_map = []
    for node in G:
        if G.nodes[node]['alive'] == [0]:
            colour_map.append('gray')
        elif G.nodes[node]['state'][t] == 0:
            colour_map.append('green')
        elif G.nodes[node]['state'][t] == 1:
            colour_map.append('red')
        elif G.nodes[node]['state'][t] == 2:
            colour_map.append('orange')
    return colour_map

def compare(G, n, lst):

    spread(G, n, lst)

    plt.plot(data(G)[0], 'g-', label='susceptible')  # plotting time evolution of number of susceptible
    plt.plot(data(G)[1], 'r-', label='active TBD')  # plotting time evolution of number of TBD
    plt.plot(data(G)[2], 'y-', label='LTB infected')  # plotting time evolution of number of high risk
    plt.xlabel('time (months)')
    plt.ylabel('population')
    plt.title('Graph showing the time evolution of the TB states of the simulated population')
    plt.legend()
    plt.show()

    nx.draw(G, node_color=colour(G, t), node_size=30)  # --> draws the graph G at end of trial
    plt.show()

if __name__ == "__main__":
    n = 100  # number of nodes at time t = 0.
    m = 3
    G = setup(n, m)
    nx.draw(G, node_color=colour(G, 0), node_size=30)  # --> draws the graph G at start of trial
    plt.show()

    none = G
    WA = G
    WB = G
    WC = G
    CA = G
    CB = G
    CC = G

    compare(none, n, [0, 'none', 0])
    #compare(WA, n, [0, 'wide', 'A'])
    #compare(WB, n, [0, 'wide', 'B'])
    #compare(WC, n, [0, 'wide', 'c'])

    #compare(CA, n, [0, 'cluster', 'A'])
    #compare(CB, n, [0, 'cluster', 'B'])
    #compare(CC, n, [0, 'cluster', 'C'])
