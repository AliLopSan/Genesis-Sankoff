import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
from treebased import TB_Network
from tools import print_cs_matrix

from mne_connectivity.viz import plot_connectivity_circle

def plot_fa_network(fas,name):
    nodes = set()
    for fa in fas:
        for node in fas[fa]:
            nodes.add(node)
    names = [n.label for n in list(nodes)]
    non_zero = []
    for char in fas:
        fa_list = list(fas[char])
        for node1 in fa_list:
            i = names.index(node1.label)
            for node2 in fa_list:
                if node1 != node2:
                    j = names.index(node2.label)
                    non_zero.append((i,j))
    con_matrix = np.zeros((len(nodes),len(nodes)))
    for i,j in non_zero:
        con_matrix[i,j] = con_matrix[i,j] + 100

    fig,axes = plot_connectivity_circle(con_matrix, names)
    fig.savefig(name,format='pdf',bbox_inches='tight')
    plt.show()
    plt.close()

f1 = open("./real_data/species_tree_ultra.pkl","rb")
S  = pkl.load(f1)
f1.close()

f2 = open("./data/interphylum_matrix.pkl","rb")
m = pkl.load(f2)
f2.close()

#Create a TB Network
N = TB_Network(S.root)
N.init_base_from_tralda(S,4)

for leaf in N.leaves():
    leaf.chars = m[leaf.label]

#Generate a basic Fitch Labeling
N.fitch_labeling()
fas = N.get_fas_by_state_change()
          
plot_fa_network(fas,"potentially_galled_characters.pdf")
