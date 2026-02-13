import itertools
from operator import mul
import sage.graphs.graph

load("lib.sage")

class Node:
    def __init__(self, parent=None, data=None, children=[]):
        self.parent = parent
        self.data = data
        self.children = children

    def add_child(self,child,i=None):
        child.parent = self
        if i != None:
            self.children.insert(child,i)
        else:
            self.children.append(child)
        return self

    def to_lst(self):
        lsts = [x.to_lst() for x in self.children]
        if self.data != None and lsts != []:
            return [self.data] +  lsts
        elif self.data != None:
            return self.data
        else:
            return lsts

    def flatten(self):
        flats = [x.flatten() for x in self.children]

        if self.data != None and flats != []:
            return [self.data] +  list(itertools.chain.from_iterable(flats))
        elif self.data != None:
            return [self.data]
        else:
            return list(itertools.chain.from_iterable(flats))

    def from_list_only_leafs(lst_or_dat):
        create = Node(None,None,[])
        if isinstance(lst_or_dat,list):
            nodes = [Node.from_list_only_leafs(x) for x in lst_or_dat]
            for x in nodes:
                x.parent = create
            create.children = nodes
        else:
            #print(lst_or_dat)
            create.data = lst_or_dat
        return create

    def from_data(dat):
        create = Node(None,dat,[])
        return create

    def left_leaf(self):
        if self.children != []:
            return self.children[0].left_leaf()
        else:
            return self

    def right_leaf(self):
        if self.children != []:
            return self.children[-1].right_leaf()
        else:
            return self

    def all_leafs(self):
        if self.children != []:
            return list(itertools.chain.from_iterable([x.all_leafs() for x in self.children]))
        else:
            return [self]

    def deep_copy(self,new_parent=None):
        create = Node(None,None,[])
        create.children = [x.deep_copy(create) for x in self.children]
        create.data = self.data
        if self.parent != None:
            assert new_parent != None
            create.parent = new_parent
        else:
            create.parent = None
        return create


def irr_on_pants(a,b,c):
    return a**2 + b**2 + c**2 != 4

def no_two_but_punct(pants):
    punct = pants["build"] == None
    sphere = pants["tr"] == None
    check_me = punct or sphere or abs(pants["tr"]) != 2
    check_inner = punct or all([no_two_but_punct(x) for x in pants["build"]])

    return check_me and check_inner

def all_inf_order(pants):
    Id = matrix([[1,0],[0,1]])
    check_me = not pants["mat"] or pants["mat"]**12 != Id
    check_inner = pants["build"] == None or all(all_inf_order(x) for x in pants["build"])
    return check_me and check_inner

def bd_trace_pants(pants):
    if pants["build"] == None:
        return None
    x = pants["build"][0]["tr"]
    y = pants["build"][1]["tr"]
    if pants["tr"]:
        z = pants["tr"]
    else:
        z = pants["build"][2]["tr"]
    return (x,y,z)

def irr_on_comps(pants):
    traces = bd_trace_pants(pants)
    check_me = pants["build"] == None or irr_on_pants(*traces)
    check_inner = pants["build"] == None or all(irr_on_comps(x) for x in pants["build"])
    return check_me and check_inner

def flatten(container):
    for i in container:
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i

def build_from_pants(pants_idx,rho,sphere=True):
    is_punct = None
    inner = not(sphere)
    if not isinstance(pants_idx,list):
        assert inner
        return {
            "build": is_punct,
            "mat": rho[pants_idx],
            "tr": rho[pants_idx].trace(),
        }
    if not(len(pants_idx) == 3 or (len(pants_idx) == 2 and inner)):
        #print(pants_idx)
        assert False,"invalid pants"
    inner_build = [build_from_pants(x,rho,False) for x in pants_idx]
    if not inner:
        flat = list(flatten(pants_idx))
        assert flat == list(range(0,max(flat)+1)), "invalid pants"
        assert inner_build[0]["mat"] * inner_build[1]["mat"] == inner_build[2]["mat"]**(-1), "invalid pants"

        return {
            "build": inner_build,
            "mat": None,
            "tr": None,
        }
    else:
        mat = inner_build[0]["mat"]*inner_build[1]["mat"]
        return {
            "build": inner_build,
            "mat": mat,
            "tr": mat.trace()
        }

def construct_binary_nested(L):
    if len(L) == 1:
        return [L[0]]

    possible_structures = []

    for i in range(1, len(L)):
        left_part = L[:i]
        right_part = L[i:]

        left_structures = construct_binary_nested(left_part)
        right_structures = construct_binary_nested(right_part)

        for l_struct in left_structures:
            for r_struct in right_structures:
                possible_structures.append([l_struct, r_struct])
    return possible_structures

def construct_ternary_root_binary_subs(L,bd=None):
    n = len(L)
    if n < 3:
        return []

    valid_roots = []

    for i in range(1, n - 1):
        for j in range(i + 1, n if not bd else i+bd+1):
            part1 = L[:i]
            part2 = L[i:j]
            part3 = L[j:]

            structs1 = construct_binary_nested(part1)
            structs2 = construct_binary_nested(part2)
            structs3 = construct_binary_nested(part3)

            for s1 in structs1:
                for s2 in structs2:
                    for s3 in structs3:
                        valid_roots.append([s1, s2, s3])

    return valid_roots

def litt_coccia(pants_idx,rho):
    pants = build_from_pants(pants_idx,rho)
    return all_inf_order(pants) and irr_on_comps(pants) and no_two_but_punct(pants)

def litt_coccia_any(rho):
    lst = range(0,len(rho))
    possible_pants = construct_ternary_root_binary_subs(lst,bd=1)
    #print(possible_pants)
    return any([litt_coccia(x,rho) for x in possible_pants])


# In[25]:



def build_special_rep(d):
    n = 12*d
    rho = rho_std(12*d)
    Br = BraidGroup(12*d)
    Br_gens = Br.gens()
    braid1 = Br_gens[1]*Br_gens[0]**2 * Br_gens[3]
    move_guys1 = Br_gens[3]**(-1) * Br_gens[4]**(-1) * Br_gens[5]**(-1)
    move_guys2 = Br_gens[5]**(-1) * Br_gens[6]**(-1)
    rep = twist_rho_mat(move_guys2*move_guys1* Br_gens[2]**3*braid1,rho=rho)

    braid_fix = reduce(mul,[Br_gens[i] for i in range(12,12*d,6)],Br_gens[0]**0)
    rep = twist_rho_mat(braid_fix,rho=rep)
    return rep

def add_punct(pant_tree):
    curr_right = pant_tree.right_leaf()
    curr_middle = pant_tree.children[1]
    new_middle = Node(pant_tree,None,[curr_middle,curr_right])
    new_right = Node(pant_tree,curr_right.data+1,[])
    pant_tree.children[1] = new_middle
    pant_tree.children[-1] = new_right

def example_pants(d):
    pants = [[0,[1,[2,3]]],[[[[[[4,5],6],7],8],9],10],11]
    new_pants = Node.from_list_only_leafs(pants)
    for i in range(0,12*(d-1)):
        add_punct(new_pants)
    return new_pants

for d in range(1,10):
    if litt_coccia(example_pants(d).to_lst(),build_special_rep(d)):
        print("Confirmed Litt-Coccia conditions for n = 12*",d)

class LabelledVertex(str):
    def __new__(cls, name, unique_id=None):
        return super(LabelledVertex, cls).__new__(cls, name)

    def __init__(self, name, unique_id=None):
        self.unique_id = unique_id if unique_id is not None else name

    def __repr__(self):
        return str(self)

def create_tree_graph(d,printer=False):
    rep = build_special_rep(d)
    assert reduce(mul,rep,matrix([[1,0],[0,1]])) == matrix([[1,0],[0,1]])
    the_pants_tree = build_from_pants(example_pants(d).to_lst(),rep)
    G = Graph(0)
    id_ctr = 0
    def _process(nod,G,parent=None):
        nonlocal id_ctr
        id_ctr += 1

        t = bd_trace_pants(nod)
        if t != None:
            print(t)
        me = (str(t) if t != None else "",id_ctr)
        G.add_vertex(me)
        if parent:
            G.add_edge(me,parent)
        if nod["build"]:
            for x in nod["build"]:
                _process(x,G,me)
    _process(the_pants_tree,G)
    return G

G = create_tree_graph(2,printer=True)
G_plot = G.plot(vertex_labels={v: v[0] for v in G.vertices()},layout='tree',tree_root=G.vertices()[0])
G_plot.save('tree_graph.png',figsize=[5*6.4,5*4.8])

