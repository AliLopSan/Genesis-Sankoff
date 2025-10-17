# version 1.2
#-----------------------------------------------------------
# eliminates the need for a labels_to_n dict and tralda_tb_dict
# when starting from a tralda tree.
#............................................................

from tralda.datastructures.Tree import TreeNode,LCA
from collections import Counter

import numpy as np

class TB_Node(TreeNode):
    def __init__(self,number_characters=1):
        super().__init__()
        self.ntype = None
        self.label = None
        self.chars = [None] * number_characters
        self.tstamp = None
        self.tb_to_tralda = None
        self.tralda_to_tb = None

class TB_Network():
    def __init__(self,root_node):
        self.root = root_node
        #Base Tree
        self.S = None
        self.transfers = [] #transfer edge list
        self.transfer_weight = dict() # transfer weight dict (key:position in transfer list)

    def preorder(self):
        def _preorder(node):
            yield node
            for child in node.children:
                yield from _preorder(child)
        if self.root:
            yield from _preorder(self.root)
        else:
            yield from []
            
    def traverse_subtree(self,u):
        yield u
        for child in u.children:
            yield from self.traverse_subtree(child)

    def postorder(self):
        def _postorder(node):
            for child in node.children:
                yield from _postorder(child)
            yield node
        if self.root:
            yield from _postorder(self.root)
        else:
            yield from []

    def inner_nodes(self):
        def _inner_nodes(node):
            if node.children:
                yield node
                for child in node.children:
                    yield from _inner_nodes(child)

        if self.root:
            yield from _inner_nodes(self.root)
        else:
            yield from []

    def leaves(self):
        def _leaves(node):
            if not node.children:
                yield node
            else:
                for child in node.children:
                    yield from _leaves(child)

        if self.root:
            yield from _leaves(self.root)
        else:
            yield from []
            
    def transfer_nodes(self):
        def _trans_nodes(node):
            if node.ntype == 'transfer':
                yield node
            else:
                for child in node.children:
                    yield from _trans_nodes(child)

        if self.root:
            yield from _trans_nodes(self.root)
        else:
            yield from []

        
            
    def tree_edges(self):
        def _edges(node):
            for child in node.children:
                yield (node, child)
                yield from _edges(child)
        
        if self.root:
            yield from _edges(self.root)
        else:
            yield from []

    def is_ancestor(self,u,v):
        current = v
        ans = False
        while current!= None:
            if current.parent == u:
                ans = True
                current = None
            else:
                current = current.parent
        return ans

    def delete_and_reconnect(self,node):
        parent = node.parent
        if not parent:
            return False
        else:
            children = [child for child in node.children]

            for child in children:
                parent.add_child_right_of(child,node)

            parent.remove_child(node)
            node.children.clear()

        return parent

    #Starts a base tree from an asymmetree species tree 
    def init_base_from_tralda(self,S,n_characters):
        from_tralda_to_tb = dict()
        for n in S.postorder():
            new = TB_Node(n_characters)
            new.ntype = 'tree'
            new.label = n.label
            new.chars = [False]*n_characters
            new.tstamp = n.tstamp
            from_tralda_to_tb[n] = new
    
            if not n.is_leaf():
                for child in n.children:
                    new.add_child(from_tralda_to_tb[child])
        self.S = TB_Network(from_tralda_to_tb[S.root])
        self.root = from_tralda_to_tb[S.root]

        self.tralda_to_tb = from_tralda_to_tb
        labels_to_n = dict()
        for v in S.postorder():
            labels_to_n[v.label]  = v
        self.tb_to_tralda = labels_to_n
        

    #Generates a labeling for the base tree
    def fitch_labeling(self):
        for node in self.postorder():
            if not node.is_leaf():
                temp = [True]*len(self.root.chars)
                for child in node.children:
                    temp = [x and y for x,y in zip(temp,child.chars)]
                node.chars = temp
                

    def sankoff_labeling(self,loss_cost,fa_cost):
        def _inner_label_cost(current_label,v,dp_table):
            label_sum = 0
            for child in v.children:
                label_cost = np.copy(dp_table[child])
                for label in range(0,2):
                    label_cost[label] = label_cost[label] + M[current_label,label]
                label_sum = label_sum + label_cost.min()
            return label_sum
        #cost matrix
        M = np.array([[0,fa_cost],[loss_cost,0]])

        for char in range(0,len(self.root.chars)):
            #Dynamic programming table initialization
            dp_table = dict()

            #labels for backtrack stage
            labels = dict()
    
            for v in self.postorder():
                dp_table[v] = np.array([0,0])

            #Sankoff's algorithm iterative version
            for v in self.postorder():
                #Leaf initialization
                if v.is_leaf():
                    if v.chars[char] == False:
                        dp_table[v][0] = 0.0
                        dp_table[v][1] = 1000000000 #not inf but something very big
                    else:
                        dp_table[v][0] = 1000000000 # same
                        dp_table[v][1] = 0.0
                else:
                    dp_table[v][0] = _inner_label_cost(0,v,dp_table)
                    dp_table[v][1] = _inner_label_cost(1,v,dp_table)
                #Storing backtrack values
                label_choice = 0
                for i in range(0,len(dp_table[v])):
                    if dp_table[v][i] < dp_table[v][label_choice]:
                        label_choice = i
                labels[v] = label_choice
                
            #Backtrack stage
            for v in self.postorder():
                if labels[v]:
                    v.chars[char] = True
                else:
                    v.chars[char] = False


    def get_types_of_transfers(self,S,n_characters):
        def _replace_w_clade_test(self,root,n_chars):
            leaves = []
            for u in self.traverse_subtree(root):
                if u.is_leaf():
                    leaves.append(u)
            temp = [True]*n_chars
            for u in leaves:
                temp = [x and y for x,y in zip(temp,u.chars)]

            if any(temp):
                return True
            else:
                return False
        clade_replaceable = []
        triangles = []
        lca_T = LCA(S)
        for donor,recipient in self.transfers:
            if donor[0] == recipient [0]:
                triangles.append((donor,recipient))

            lca = lca_T(self.tb_to_tralda[donor[1]],self.tb_to_tralda[recipient[1]])
            if _replace_w_clade_test(self,self.tralda_to_tb[lca],n_characters):
                clade_replaceable.append((donor,recipient))

        return(clade_replaceable,triangles)

    def get_fas_by_state_change(self):
        first_apps = dict()

        #Initialize dict
        for i in range(0,len(self.root.chars)):
            first_apps[i] = []

        #Change counter
        for node in self.postorder():
            if node.parent!= None:
                for i in range(0,len(node.chars)):
                    if node.chars[i] == True and node.parent.chars[i] == False:
                        temp = list(first_apps[i])
                        temp.append(node)
                        first_apps[i] = temp
        return(first_apps)

    def greedy_completion(self,fa_list,S):
        def _sort_fas(fa_list):
            return(sorted(fa_list, key = lambda x:x.tstamp,reverse=True))
        def _create_transfer(lca_T,underlying_outt,donor,recipient):
            self.transfers.append((donor,recipient))
            self.transfer_weight[len(self.transfers) - 1] = 0
            lca = lca_T(self.tb_to_tralda[donor[1].label],self.tb_to_tralda[recipient[1].label])
            current = donor[1].parent
            while current!= self.tralda_to_tb[lca]:
                underlying_outt[current].add(recipient[1])
                current = current.parent
        def _look_for_donor_edge(donor,recipient):
            if donor[1].is_leaf():
                return donor
            else:
                for child in self.traverse_subtree(donor[1]):
                    if child.parent.tstamp > recipient[1].tstamp and child.tstamp <= recipient[1].tstamp:
                        return ((child.parent,child))
        #Auxiliary LCA structure
        lca_T = LCA(S)

        #Start a dictionary of underlying outtransfers for each node
        underlying_outt = dict()
        for node in self.postorder():
            underlying_outt[node] = set()
            
        for character in range(0,len(self.root.chars)):
            if len(fa_list[character]) > 1:
                X = _sort_fas(fa_list[character])
                for i in range(0,len(X)-1):
                    donor =(X[i].parent,X[i])
                    recipient = (X[i+1].parent,X[i+1])
                    if (donor,recipient) in self.transfers:
                        pos = self.transfers.index((donor,recipient))
                        self.transfer_weight[pos] =+1
                    else:
                        if X[i+1] not in underlying_outt[X[i]]:
                            donor = _look_for_donor_edge(donor,recipient)
                            _create_transfer(lca_T,underlying_outt,donor,recipient)

                            
    #After the greedy completion, output the transfer info
    def print_transfer_hw_info(self):
        if len(self.transfers) == 0:
            print("\t There are no transfer arcs in this network :) ")
        else:
            normal   = 0
            h_list = []
            donors_list = []
            recipients_list = []
            max_trans = -1
            for i in range(0,len(self.transfers)):
                donor,recipient = self.transfers[i]
                donors_list.append(donor)
                recipients_list.append(recipient)
                if self.transfer_weight[i] > max_trans:
                    max_trans = self.transfer_weight[i]
                if self.transfer_weight[i] > 0:
                    h_list.append((self.transfers[i],self.transfer_weight[i]))
                else:
                    normal +=1

            print("\t There are ",len(h_list)," transfer highways vs ",normal," normal transfers")
            print("\t Maximum weight on highway is ", max_trans)
            occ_count_don = Counter(donors_list)
            occ_count_rec = Counter(recipients_list)
            print("\t The most common donor edge: ",occ_count_don.most_common(1)[0])
            print("\t The most common recipient is:",occ_count_rec.most_common(1)[0])
            
            
    def genesis_labeling(self, loss_cost, fa_cost):
        def _base_case_value(dp_table, v, label, gain, char):
            inf = 1000000000  # not inf but something very big
            if v.chars[char] == False:
                if label == 0 and gain == 0:
                    return 0
                else:
                    return inf
            else:
                if label == 1:
                    return 0
                else:
                    return inf

        def _inner_case_value(dp_table, v, current_label, current_gain, M):
            # Special case: V[v,0,1] - Loss state (0) with gain origin (1)
            # This case fixes one child as origin and leaves the rest to fate
            if current_label == 0 and current_gain == 1:
                child_sum = []
                for child in v.children:
                    origin_penalty = []
                    for l in [0, 1]:
                        origin_penalty.append(dp_table[child][l, 1])

                    brother_sum = 0
                    for brother in v.children:
                        if brother != child:
                            label_cost = []
                            for l in [0, 1]:
                                # Brothers are non-origin, so they are not allowed to be gain origins (gain=0)
                                # The transition cost is M[0, l] since parent is labeled 0 (loss).
                                label_cost.append(dp_table[brother][l, 0] + M[0, l])
                            brother_sum = brother_sum + min(label_cost)
                    child_sum.append(min(origin_penalty) + brother_sum)
                return min(child_sum)
            else:
                # All other states V[v, 1, 0], V[v, 1, 1], V[v, 0, 0]
                # These states do not require a special 'origin' child, so all children contribute
                cost = 0
                for child in v.children:
                    label_costs = []
                    # Transition from label l at child to current_label at v. Child must be non-origin (gain=0)
                    # The original function had l in [1, 0] which is reversed from standard practice, 
                    # but kept here for consistency with the provided snippet.
                    for l in [1, 0]:
                        label_costs.append(dp_table[child][l, 0] + M[current_label, l])
                    cost = cost + min(label_costs)
                return cost

        # Backtracking auxiliary function for choosing the best origin child when V[v, 0, 1] is chosen
        def _find_best_child_origin(dp_table, v, M):
            argmin = []
            children = list(v.children)
            # M[0, l] is used because the parent label is 0 (loss) in this V[v, 0, 1] case.
            for child in v.children:
                origin_penalty = []
                for l in [0, 1]:
                    origin_penalty.append(dp_table[child][l, 1] + M[0, l])
                
                brother_sum = 0
                for brother in v.children:
                    if brother != child:
                        label_cost = []
                        for l in [0, 1]:
                            label_cost.append(dp_table[brother][l, 0] + M[0, l])
                        brother_sum = brother_sum + min(label_cost)
                
                argmin.append(min(origin_penalty) + brother_sum)
            
            index_min = np.argmin(argmin)
            return children[index_min]

        # cost matrix: M[parent_label, child_label]
        # label 0: False (absence), label 1: True (presence/gain)
        # M[0, 1]: transition 0 -> 1 (gain) cost is fa_cost
        # M[1, 0]: transition 1 -> 0 (loss) cost is loss_cost
        M = np.array([[0.0, fa_cost], [loss_cost, 0.0]])
        dp_table = dict()

        # Main algorithm
        for char in range(0, len(self.root.chars)):
            dp_table = dict()
            
            # 1. Forward (Postorder) Pass: Fill DP table
            for v in self.postorder():
                dp_table[v] = dict()
                # label 0: False (absence/loss), label 1: True (presence/gain)
                # gain 0: no gain origin in subtree, gain 1: gain origin in subtree
                for label in [0, 1]:
                    for gain in [0, 1]:
                        if v.is_leaf():
                            dp_table[v][label, gain] = _base_case_value(dp_table, v, label, gain, char)
                        else:
                            dp_table[v][label, gain] = _inner_case_value(dp_table, v, label, gain, M)

            # 2. Backtracking (Preorder) Pass: Determine labels and origin
            origin = None # Tracks the node that is the 'origin' of the gain for the current character
            
            for v in self.preorder():
                parent_char_is_present = v.parent and v.parent.chars[char]
                
                if v == self.root:
                    # Root: Must be a gain origin (gain=1) or not present (label=0, gain=0)
                    # The root's state is chosen between V[root, 1, 1] and V[root, 0, 1]
                    # V[root, 0, 1] allows the gain to originate *below* the root.
                    if dp_table[v][1, 1] <= dp_table[v][0, 1]:
                        v.chars[char] = True # V[root, 1, 1] chosen
                        origin = self.root # Root is the origin
                    else:
                        v.chars[char] = False # V[root, 0, 1] chosen
                        # choose the origin child
                        origin = _find_best_child_origin(dp_table, v, M)
                
                elif v == origin:
                    # Current node is the gain origin. Must choose between V[v, 1, 1] and V[v, 0, 1]
                    # If V[v, 0, 1] is chosen, the true origin must be a child of v.
                    if dp_table[v][1, 1] <= dp_table[v][0, 1]:
                        v.chars[char] = True
                    else:
                        v.chars[char] = False
                        # choose the origin child from v's subtree
                        if not v.is_leaf():
                            origin = _find_best_child_origin(dp_table, v, M)
                
                else:
                    # Current node is not the gain origin. Must choose between V[v, 1, 0] and V[v, 0, 0]
                    # Only the V[v, label, 0] states are allowed.
                    if dp_table[v][1, 0] <= dp_table[v][0, 0]:
                        v.chars[char] = True
                    else:
                        v.chars[char] = False
            

    

        
            
    
        
        
