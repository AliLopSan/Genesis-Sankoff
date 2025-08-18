#................................................
#      GENERATE CHARACTER-STATE MATRICES
#................................................
#                  version 1.0
import requests
import random as rd
import pickle as pkl
import pandas as pd
import asymmetree.treeevolve as te
from treebased import TB_Network
from alive_progress import alive_bar


#INPUT: S species tree with tralda, N TB_Network 
#OUTPUT: N annotated with presence/absence data

#................................................
#  Generate a random state list
#...............................................
def generate_random_matrix(S,n_chars):
    character_state = dict()
    for species in S.leaves():
        temp = []
        for i in range(0,n_chars):
            coin_flip = rd.randint(0,1)
            if coin_flip == 0:
                temp.append(False)
            else:
                temp.append(True)
        character_state[species.label] = list(temp)
    
    return(character_state)
    
    
#.................................................
# Generate a matrix with simulated data
# Additional input: number of characters
# Saves the ground-truth highway graph
#................................................


def assign_states_from_gt(char,ogt,H):
    def _get_reconc_transfer_edge(node):
        if node.event != "H":
            return KeyError(f'{node} is not an HGT node')
        for v in node.children:
            if hasattr(v,'transferred') and v.transferred:
                if v.event == "S":
                    parent = H.tb_to_tralda[v.reconc].parent
                    recipient = (parent.label,v.reconc)
                else:
                    recipient = v.reconc
                yield(node.reconc,recipient)
    transfer_edges_set = set(H.transfers)
    inner = [node for node in ogt.inner_nodes()]
    root  = rd.choice(inner)

    for n in ogt.traverse_subtree(root):
        if n.event != 'H':
            node_in_tb = H.tralda_to_tb[H.tb_to_tralda[n.reconc]]
            node_in_tb.chars[char] = True
        else:
            for edge in _get_reconc_transfer_edge(n):
                if edge in transfer_edges_set:
                    position =H.transfers.index(edge)
                    H.transfer_weight[position] = H.transfer_weight[position] + 1
                else:
                    H.transfers.append(edge)
                    transfer_edges_set.add(edge)
                    H.transfer_weight[len(H.transfers) - 1] = 0


def generate_from_asymm(S,N,n_chars,params_dict,name):
    character_state = dict()
    #Generate highway graph
    H = TB_Network(S.root)
    H.init_base_from_tralda(S,n_chars)

    #Generate inner labeling
    for i in range(0,n_chars):
        #Generation of true gene tree (tgt)
        if params_dict['HGT_type'] == "exponential":
            strength = rd.uniform(0,1)
            tgt = te.dated_gene_tree(S,
                                   dupl_rate=params_dict["dup_rate"],
                                   loss_rate=params_dict["loss_rate"],
                                     hgt_rate=params_dict["hgt_rate"],
                                   additive_transfer_distance_bias='exponential',
                                   transfer_distance_bias_strength=strength,
                                     prohibit_extiction='per_familiy')
        elif params_dict['HGT_type'] == "inverse":
            strength = rd.uniform(0,1)
            tgt = te.dated_gene_tree(S,
                                   dupl_rate=params_dict["dup_rate"],
                                   loss_rate=params_dict["loss_rate"],
                                     hgt_rate=params_dict["hgt_rate"],
                                   additive_transfer_distance_bias='inverse',
                                   transfer_distance_bias_strength=strength,
                                     prohibit_extiction='per family')
        else:
            tgt = te.dated_gene_tree(S,
                                     dupl_rate=params_dict["dup_rate"],
                                     loss_rate=params_dict["loss_rate"],
                                     hgt_rate=params_dict["hgt_rate"],
                                     prohibit_extiction='per_familiy')

        #print("\t Original gene tree: ",tgt.to_newick())
        ogt = te.prune_losses(tgt)
        #print("\t Prunned tree :",ogt.to_newick())
        assign_states_from_gt(i,ogt,H)

    #Output resulting highway graph
    H_name = name + '_highway_graph.pkl'
    file = open(H_name,'wb')
    pkl.dump(H,file)
    file.close()

    for leaf in H.leaves():
        character_state[leaf.label] = leaf.chars

    return (character_state,H)
            
                
        



#.................................................
# Generate a matrix with real-life data using KEGG
# INPUT: Two lists
#    - A list of species ids of the form (T000) that
#     correspond to the species at the leaves
#    - A list of K0s of the form (K0789475) that
#     correspond to the characters
#.................................................
def generate_from_KEGG_v1(org_list,ko_list):
    char_state = dict()
    print("\t Building char-state matrix: ")
    with alive_bar(len(org_list)) as bar:
        for org in org_list:
            states = []
            for ko in ko_list:
                url = "https://rest.kegg.jp/link/" + org + "/ko:" + ko
                r = requests.get(url)
                if (r.status_code == 400):
                    print("ERROR: Bad request to server")
                    sys.exit()
                else:
                    if (len(r.text) == 1):
                        states.append(0)
                    else:
                        states.append(1)
            char_state[org] = states
            bar()
        
    return(char_state)
