import networkx as nx
from enum import Enum
from .monodromy_rep import MonodromyRep
from itertools import repeat
from multiprocessing import Pool

class OrbitGraph:
    num_gens = 0
    gens = []
    start = None
    action = None
    ident = None
    digraph = None

    def __init__(self, _gens, _start, _action = MonodromyRep.action, **args):
        self.gens = _gens
        self.action = _action
        self.ident = _gens[0] ** 0
        self.num_gens = len(_gens)

        # Maps hash → list of canonical node objects. Ensures each unique node
        # is represented by a single Python object so that NetworkX dict lookups
        # hit the fast `is`-identity path rather than calling __eq__ every time.
        self._node_index = {}

        self.digraph = nx.MultiDiGraph()

        _start = self._intern(_start)
        self.start = _start
        self.digraph.add_node(_start)

    class DIR(Enum):
        FORWARD = 1
        BACK = 2

    def _intern(self, node):
        """Return the canonical object for *node*.

        On first encounter the node itself is stored as canonical.  On
        subsequent encounters with a mathematically-equal node the previously
        stored object is returned, so NetworkX dict lookups can use the fast
        `is`-identity path instead of calling __eq__ each time.
        """
        h = hash(node)
        for existing in self._node_index.get(h, ()):
            if existing is node or existing == node:
                return existing
        self._node_index.setdefault(h, []).append(node)
        return node

    def advance_node_gen_pure(self,n,g,direction):
        if direction == self.DIR.FORWARD:
            m = self.action(g,n)
            return (m, (n,m,str(g)))
        elif direction == self.DIR.BACK:
            m = self.action(g**int(-1),n)
            return (m,(m,n,str(g)))
        else:
            assert False, "Direction Must be a DIR"

    def advance_node_gen_pure_check(self,n,g,direction):
        (node,edge) = self.advance_node_gen_pure(n,g,direction)
        if self.digraph.has_edge(*edge):
            return None
        else:
            return (node,edge)

    def advance_node_gen(self, n, g, direction):
        g_key = str(g)
        if direction == self.DIR.FORWARD:
            # Skip the twist entirely if this generator's forward edge already exists.
            if any(g_key in edges for edges in self.digraph.succ[n].values()):
                return
            m = self._intern(self.action(g, n))
            self.digraph.add_node(m)
            self.digraph.add_edge(n, m, key=g_key)
        else:
            # Skip the twist entirely if this generator's backward edge already exists.
            if any(g_key in edges for edges in self.digraph.pred[n].values()):
                return
            m = self._intern(self.action(g ** int(-1), n))
            self.digraph.add_node(m)
            self.digraph.add_edge(m, n, key=g_key)

    def advance_node_dir(self,n,direction,par=False):
        if par:
            with Pool() as pool:
                func = self.advance_node_gen_pure_check
                lst = pool.starmap(func, zip(repeat(n), self.gens, repeat(direction)))
                self.add_all_nodes_edges(lst)
        else:
            for g in self.gens:
                self.advance_node_gen(n,g,direction)

    def advance_node_for_par(self,n,direction):
        lst = [self.advance_node_gen_pure_check(n,g,direction) for g in self.gens]
        return lst

    def add_all_nodes_edges(self,lst):
        for x in lst:
            if x is not None:
                (node, e) = x
                node = self._intern(node)
                (v, w, k) = e
                # Replace whichever endpoint is the new node with its canonical form.
                if w is not v:
                    e = (v, node, k) if node == w else (node, w, k)
                self.digraph.add_node(node)
                self.digraph.add_edge(*e)

    def advance_lst_dir(self, inp_lst, direction, par=False):
        lst = [node for node, i in inp_lst if i != self.num_gens]
        if par:
            with Pool() as pool:
                lsts = pool.starmap(self.advance_node_for_par, zip(lst, repeat(direction)), chunksize=(int(len(lst)/20)+1))
                pool.close()
                pool.join()
                for to_add in lsts:
                    self.add_all_nodes_edges(to_add)
        else:
            for node in lst:
                self.advance_node_dir(node, direction)

    def advance_forward(self):
        lst = self.digraph.out_degree()
        self.advance_lst_dir(lst, self.DIR.FORWARD)

    def advance_back(self):
        lst = self.digraph.in_degree()
        self.advance_lst_dir(lst, self.DIR.BACK)

    def advance(self):
        lst1 = list(self.digraph.out_degree())
        lst2 = list(self.digraph.in_degree())

        self.advance_lst_dir(lst1, self.DIR.FORWARD)
        self.advance_lst_dir(lst2, self.DIR.BACK)

    def advance_until(self, N = -1):
        if N == -1:
            last_size = -1
            this_size = len(self.digraph.nodes())
            while last_size != this_size:
                last_size = this_size
                self.advance()
                this_size = len(self.digraph.nodes())
        else:
            for i in range(0, N):
                self.advance()

    def show(self, node_labels=False, edge_labels=True,colors={}):
        return None
