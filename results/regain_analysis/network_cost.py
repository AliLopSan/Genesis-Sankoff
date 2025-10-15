import pickle as pkl
import pandas as pd
import sys

def network_cost_by_character(network,loss_penalty,trans_penalty):
    losses_per_character = []
    fas_per_character = []
    for char in range(0,len(N.root.chars)):
        transfer_cost = 0
        loss_cost     = 0

        for u,v in N.tree_edges():
            label_u = u.chars[char]
            label_v = v.chars[char]
            if label_u != label_v:
                if label_u == True and label_v == False:
                    loss_cost+=1
                else:
                    transfer_cost+=1
            if transfer_cost > 1:
                transfer_cost = transfer_cost - 1
        loss_cost = loss_cost*loss_penalty
        transfer_cost = transfer_cost*trans_penalty
        losses_per_character.append(loss_cost)
        fas_per_character.append(transfer_cost)

    results = pd.DataFrame()
    results['losses'] = losses_per_character
    results['transfers'] = fas_per_character
    results['network_cost'] = results['losses'] + results['transfers']
    print(results.head())
    return(results)
        

file = open("sankoff_network_hw1.pkl","rb")
N = pkl.load(file)
N_cost = 0
file.close()

r = network_cost_by_character(N,1.0,1.0)
        
