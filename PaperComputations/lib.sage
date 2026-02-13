

# # Setup Cells

# ## Braid Group Action Setup

### FREE GROUP -> SL2Z and Rep Defns ### 
import math

def free_to_matrix(f, rho):
    res = rho[0]**0
    for i in f.Tietze():
        res = res * rho[abs(i)-1]**sgn(i)
    return res

A = matrix([[1,-1],[0,1]])
B = matrix([[1,0],[1,1]])
def rho_std(n):
    rho = []
    AA = matrix([[1,-1],[0,1]])
    BB = matrix([[1,0],[1,1]])
    for i in range(0,floor(n/2)):
        rho.append(AA)
        rho.append(BB)
    if n % 2 == 1:
        rho.append(AA)
    return rho

def conj_rho(S, rho):
    return list(map(lambda x: S * x * S**(-1), rho))

def rho_traces3(rho):
    [X,Y,Z] = rho
    return [(X*Y).trace(),(Y*Z).trace(),(X*Z).trace(),(X*Y*Z).trace()]





## Testing for Conj Reps ##

def test_conj_reps(lst_reps, lst_conjs,with_idx=False):
    lst_seen = []
    lst_with_idx = []
    helpful = []
    for (j,rho) in enumerate(lst_reps):
        conjs_of_rho = [(S,lst_seen.index(conj_rho(S,rho))) for S in lst_conjs if conj_rho(S,rho) in lst_seen]
        helpful = helpful + conjs_of_rho
        if with_idx:
            for (_,k) in conjs_of_rho:
                lst_with_idx[k][2].append(j)
        if not (rho in lst_seen) and not (conjs_of_rho):
            lst_seen.append(rho)
            if with_idx:
                lst_with_idx.append((rho,j,[]))

    if with_idx:
        return (helpful,lst_with_idx)
    else:
        return (helpful, lst_seen)

def lots_of_sl2(N):
    if N == 0:
        id_mat = matrix([[1,0],[0,1]])
        return [id_mat]

    sl2_gens = [A,B,A**(-1),B**(-1)]
    return [C*X for C in lots_of_sl2(N-1) for X in sl2_gens]

def lots_of_braid(gens,N):
    if N == 0:
        return [gens[0]**0]

    gens_with_inv = [b**e for e in [1,-1] for b in gens]
    return [C*X for C in lots_of_braid(gens,N-1) for X in gens_with_inv]





### Braid Group acting on Reps + Stabilizer Test ###

def twist_rho_mat(b,rho=None):
    if rho == None:
        rho = rho_std(b.strands())
    new_rho = rho.copy()
    for i in list(b.Tietze())[::-1]:
        j = abs(i)
        if i < 0:
            tmp1 = new_rho[j-1]
            tmp2 = new_rho[j]
            new_rho[j-1] = tmp2
            new_rho[j] = tmp2**(-1) * tmp1 * tmp2
        else:
            tmp1 = new_rho[j-1]
            tmp2 = new_rho[j]
            new_rho[j-1] = tmp1*tmp2*(tmp1**(-1))
            new_rho[j] = tmp1
    return new_rho

def twist_rho(b,rho = None):
    if rho == None:
        rho = rho_std(b.strands())
    the_gens = FreeGroup(b.strands()).gens()
    return [free_to_matrix(the_gens[i]*b,rho) for i in range(0,b.strands())]

@parallel
def braid_in_img(b, rho = None):
    if rho == None:
        rho = rho_std(b.strands())
    twisted = twist_rho(b,rho)
    return rho == twisted

def braid_in_img_p(b,p=2,rho=None):
    if rho == None:
        rho = rho_std(b.strands())
    twisted = twist_rho(b,rho)
    return [mat % p for mat in rho] == [mat % p for mat in twisted]

def braid_in_core_ex(b,reps):
    for rho in reps:
        if not braid_in_img(b,rho):
            return (False,rho)
    return (True, None)

def braid_in_core(b, reps):
    return all([braid_in_img(b,rho) for rho in reps])

def find_core_tech(reps,gens,N, not_in_core_ex = None):
    ident = gens[0] ** 0
    if not_in_core_ex == None:
        not_in_core_ex = []
    not_in_core = [b for (b,_) in not_in_core_ex]

    cands = []
    for b in lots_of_braid(gens,N):
        if b in not_in_core or b == ident:
            continue

        (_, rep) = braid_in_core_ex(b,reps)
        if rep:
            not_in_core_ex.append((b,rep))
            not_in_core.append(b)
        else:
            cands.append(b)

    return (cands, not_in_core, not_in_core_ex)

def find_core_cand(reps,gens,N):
    ident = gens[0] ** 0
    cands = []
    for b in lots_of_braid(gens,N):
        if b != ident and braid_in_core(b,reps):
            cands.append(b)
    return cands


# ## Dehn Twists in T^2 version




import itertools

def scc(v):
    x = math.gcd(v[0,0],v[1,0])
    w = 1/x * v
    if w[0,0] > 0:
        return w
    elif w[0,0] == 0:
        return matrix([[0],[1]])
    else:
        return -w

def scc_from_matrix(A,simplify=True):
    e1 = matrix([[1],[0]])
    e2 = matrix([[0],[1]])

    c1 = A*e1 - e1
    c2 = A*e2 - e2
    if c1 != 0:
        return scc(c1) if simplify else c1
    elif c2 != 0:
        return scc(c2) if simplify else c2
    else:
        return None

def matrix_from_scc(c):
    e1 = matrix([[1],[0]])
    e2 = matrix([[0],[1]])
    v1 = e1 + alg_int_number(e1,c)*c
    v2 = e2 + alg_int_number(e2,c)*c
    return matrix([v1.transpose()[0], v2.transpose()[0]]).transpose()

def matrix_from_pair(p,q):
    return matrix([[1+p*q, -p**2], [q**2, 1-p*q]])

def alg_int_number(v1,v2):
    M = matrix([[v1[0,0],v2[0,0]],[v1[1,0],v2[1,0]]])
    return M.det()

def int_number(v1,v2):
    return abs(alg_int_number(v1,v2))

def scc_rep(rho,simplify=True):
    return [scc_from_matrix(X,simplify=simplify) for X in rho]

A = matrix([[1,-1],[0,1]])
B =  matrix([[1,0],[1,1]])
def adj_int_numbers(rho):
    scc_rho = scc_rep(rho)
    return [int_number(v1,v2) for (v1,v2) in itertools.pairwise(scc_rho)]

def all_int_numbers(rho):
    scc_rho = scc_rep(rho)
    return [int_number(v1,v2) for (v1, v2) in itertools.combinations(scc_rho,2)]





H = BraidGroup(2)
H_gens = H.gens()
var('y','z')

inf = Infinity
def rat_to_curve(p):
    if p == Infinity:
        return matrix([[1],[0]])
    else:
        a = p.numerator()
        b = p.denominator()
        return matrix([[a], [b]])
def curve_to_rat(c):
    if c[1][0] == 0:
        return Infinity
    else:
        return c[0][0]/c[1][0]

def apply_braid_rat(b,ps):
    rep = [matrix_from_scc(rat_to_curve(p)) for p in ps]
    new_rep = twist_rho(b,rho=rep)
    sccs_rep = scc_rep(new_rep)
    return [curve_to_rat(c) for c in sccs_rep]


apply_braid_rat(H_gens[0]**(-1), [1,Infinity])


# ## Orbit Graph Setup




# Works when rho,other are given by dehn twists and generate SL2Z
def check_conj(rho,other,n=None,make_scc=scc_rep,printer=False):
    if n == None:
        n = len(rho)
    rho_scc = make_scc(rho)
    other_scc = make_scc(other)
    fst = rho_scc[0]
    def _dim_check(v1,v2):
        return span(v1.columns() + v2.columns()).dimension() == 2

    cand_lst = [(rho_scc[i],i) for i in range(1,n) if _dim_check(fst,rho_scc[i])]
    assert cand_lst != [], "rho_scc does not span"
    snd,j = cand_lst[0]
    fst_oth = other_scc[0]
    snd_oth = other_scc[j]
    if printer:
        print(fst)
        print(snd)
        print(fst_oth)
        print(snd_oth)
    if not _dim_check(fst_oth,snd_oth):
        return False,None
    def _with_these(_fst,_fst_oth,_snd,_snd_oth):
        C1 =  matrix([[_fst[0,0], _snd[0,0]], [_fst[1,0],_snd[1,0]]])
        C2 =  matrix([[_fst_oth[0,0], _snd_oth[0,0]], [_fst_oth[1,0],_snd_oth[1,0]]])
        cand_conj = C2 * C1**(-1)
        if printer:
            print(cand_conj)

        if cand_conj.determinant() != 1:
            cand_conj[:,0] = -cand_conj[:,0]
        if other == [cand_conj * m * cand_conj**(-1) for m in rho]:
            return True,cand_conj
        else:
            return False,None
    sgns = [(1,1),(-1,1),(1,-1),(-1,-1)]
    blegh = [_with_these(fst,e1*fst_oth,snd,e2*snd_oth) for e1,e2 in sgns]
    blub = [b for b in blegh if b[0]]
    if blub == []:
        return False,None
    else:
        return blub[0]






def construct_orbit_hood(m, gens = None, rho = None,N=1000,M = 20, breadth=True, one_index_gen_names=False):
    def _stupid_conj(x,y):
        if x == y:
            return True,matrix([[1,0],[0,1]])
        else:
            return False,None
    return construct_orbit_hood_conj(m,gens,rho=rho,N=N,M=M,breadth=breadth,one_index_gen_names=one_index_gen_names,conj_checker=_stupid_conj)

def construct_orbit_hood_conj(m, gens = None, N=1000,M = 20, breadth=True, one_index_gen_names=False,rho=None,conj_checker=check_conj):
    if rho == None:
        rho = rho_std(m)
    if gens == None or gens==[]:
        if one_index_gen_names:
            BGroup = BraidGroup(m, ['σ' + str(i) for i in range(1,m)])
        else:
            BGroup = BraidGroup(m, 'σ')
        gens = BGroup.gens()
    else:
        BGroup = gens[0].parent()
    Gr = DiGraph(loops=True,multiedges=True)

    seen = [rho]
    new_vtxs = [(rho,0)]
    Gr.add_vertex()
    num_loops = 0
    i = 1
    completed_verts = []
    while new_vtxs and (i < N and num_loops < M):
        if breadth:
            (v,k) = new_vtxs.pop(0)
        else:
            (v,k) = new_vtxs.pop()
        for g in gens:
            w = twist_rho(g, v)
            w_seen = [ell for ell,r in enumerate(seen) if conj_checker(w,r)[0]]
            if w_seen == []:
                seen.append(w)
                new_vtxs.append((w,i))
                j = i
                Gr.add_vertex()
                i += 1
            else:
                j = w_seen[0]
            if not(Gr.has_edge(k,j,str(g))):
                Gr.add_edge(k,j, label=str(g))
            w = twist_rho(g**(-1), v)
            w_seen = [ell for ell,r in enumerate(seen) if conj_checker(w,r)[0]]
            if w_seen == []:
                seen.append(w)
                new_vtxs.append((w,i))
                j = i
                Gr.add_vertex()
                i += 1
            else:
                j = w_seen[0]
            if not(Gr.has_edge(j,k,str(g))):
                Gr.add_edge(j,k, label=str(g))
        completed_verts.append(k)
        num_loops += 1
    return (Gr, gens, seen, i >= N or num_loops >= M,completed_verts,BGroup)

def display_orbit(graph, gens, completed_verts = None, loopless=True, **kwargs):#, save_pos=False,layout='spring',iters=5,heights={}):
    m = len(gens)
    vert_dict = None
    if completed_verts:
        vert_dict = {}
        vert_dict[(0,1,0)] = completed_verts
        vert_dict[(1,0,0)] = [v for v in graph.vertices() if v not in completed_verts]
    bow = rainbow(m)
    color_dict = {str(gens[i]):bow[i] for i in range(0,m) }
    return graph.plot(edge_colors = graph._color_by_label(color_dict),
                      vertex_colors = vert_dict,
                      **kwargs)
                      #layout=layout,save_pos=save_pos, iterations=iters,heights=heights)


def braid_from_path(path,gens):
    lst_of_elts = [gens[int(z[1:])] for (x,y,z) in path]
    res = gens[0] ** (0)
    for x in lst_of_elts:
        res = x * res
    return res


def from_Tietze_to_genlst(g,gens):
    tietze = g.Tietze()
    return [(gens[abs(i)-1], math.copysign(1,i)) for i in tietze]

def from_Tietze_to_gens(g, gens):
    return [h**j for (h,j) in from_Tietze_to_genlst(g,gens)] 

def from_Tietze(lst,gens):
    res = gens[0]**0
    for i in lst:
        j = 1 if i > 0 else -1
        res = res * (gens[abs(i)-1] ** j)
    return res

def edge_from_gen(Gr, v, gen, k=1):
    if k == 1:
        cands = [(i,j,e) for (i,j,e) in Gr.outgoing_edges(v) if e == str(gen)]
    elif k == -1:
        cands = [(i,j,e) for (i,j,e) in Gr.incoming_edges(v) if e == str(gen)]

    if cands:
        return cands[0]
    else:
        if k == 1:
            print(Gr.outgoing_edges(v),"\n")
        else:
            print(Gr.incoming_edges(v),"\n")
        print(v,"\n",gen,"\n",k,"\n",cands,"\n\n")
        return None

def path_from_genlst_noswp(Gr, lst, v):
    w = v
    result = []
    for (gen,i) in lst:
        e = edge_from_gen(Gr,w,gen,i)
        result.append((e,i))
        (new_w_back, new_w_forward, _) = e
        w = new_w_forward if i == 1 else new_w_back
    return result

def path_from_genlst(Gr, lst, v):
    def swap_if_neg(t, j):
        if j == 1:
            return t
        if j == -1:
            return (t[1], t[0], t[2] + "^(-1)")
    return [swap_if_neg(e,i) for (e,i) in path_from_genlst_noswp(Gr,lst,v)]

def path_from_braid(Gr, gens, v, g):
    genlst = from_Tietze_to_genlst(g,gens)
    genlst.reverse()
    return path_from_genlst(Gr, genlst, v)

def antinodal_at(Gr, v, gen):
    (i, j, e) = edge_from_gen(Gr,v,gen)
    if j == v:
        return None
    [(e1,_),(e2,_),(e3,_)] = path_from_genlst_noswp(Gr, 3*[(gen,1)],v)
    (k,l,ee) = e3
    if v == l:
        return [e1,e2,e3]
    return None

def equinodal_at(Gr, v, gen):
    (i,j,e) = edge_from_gen(Gr, v, gen)
    if j == v:
        return (i,j,e)
    else:
        return None

def path_conj_braid(Gr, gens, v, w):
    a_path = Gr.all_paths(v,w,report_edges=True,labels=True)[0]
    return braid_from_path(a_path,gens)

def get_equi_antinodals(Gr,v,gens):
    equis = []
    antis = []
    ident = gens[0]**0

    for w in Gr.vertices():
        b = ident if v == w else path_conj_braid(Gr, gens, v, w)
        for gen in gens:
            if antinodal_at(Gr,w,gen):
                antis.append(b**(-1)* gen**3 *b)
            elif equinodal_at(Gr,w,gen):
                equis.append(b**(-1)*gen*b)
    return (equis, antis)



# ## Huge Orbit Setup




def twist_all(g,reps):
    return [twist_rho(g,v) for v in reps]

def construct_huge_orbit(m, reps, gens = None,N=1000,M = 20, breadth=True):
    if gens == None:
        BGroup = BraidGroup(m, 'σ')
        gens = BGroup.gens()
    Gr = DiGraph(loops=True)

    seen = [reps]
    new_vtxs = [(reps,0)]
    Gr.add_vertex()
    num_loops = 0
    i = 1
    completed_verts = []
    while new_vtxs and (i < N and num_loops < M):
        if breadth:
            (v,k) = new_vtxs.pop(0)
        else:
            (v,k) = new_vtxs.pop()
        for g in gens:
            w = twist_all(g, v)
            if not (w in seen):
                seen.append(w)
                new_vtxs.append((w,i))
                j = i
                Gr.add_vertex()
                i += 1
            else:
                j = seen.index(w)

            Gr.add_edge(k,j, label=str(g))
        completed_verts.append(k)
        num_loops += 1
    return (Gr, gens, seen, i >= N or num_loops >= M,completed_verts)


# ## Finding Generators / Misc Group Stuff

# In[10]:


import itertools as it

### Finite Generation from Finite Index ###

def same_left_coset(a,b, pred = braid_in_img):
    return pred(b**(-1) * a)
def same_right_coset(a,b, pred = braid_in_img):
    return pred(a * (b**(-1)))

def find_in_rcosets(g, rcosets, pred=braid_in_img):
    #print("    Finding for " + str(g))
    try:
        guy = next(filter(lambda x: same_right_coset(g,x,pred), rcosets))
    except:
        raise ValueError("Could not find a coset for " + str(g))
    #print("    Found: " + str(guy))
    return guy

def gens_one_rcoset(gens, cos_rep, rcosets,pred=braid_in_img):
    # print("### Working on Coset " + str(cos_rep) + " ###")
    return [cos_rep * g * (find_in_rcosets(cos_rep*g,rcosets,pred))**(-1) for g in gens]

def gens_from_rcosets(gens, rcosets, pred=braid_in_img):
    return [x for cos_rep in rcosets for x in gens_one_rcoset(gens,cos_rep,rcosets,pred)]

def rcosets_from_lcosets(lcosets):
    return list(map(lambda x: x**(-1), lcosets))

def gens_from_lcosets(gens,lcosets,pred=braid_in_img):
    rcosets = rcosets_from_lcosets(lcosets)
    return gens_from_rcosets(gens,rcosets,pred)

def build_lcosets(Gr, gens):
    lcosets = [gens[0] * gens[0]**(-1)]
    for v in Gr.vertices():
        if v != 0:
            path = next(Gr.all_paths_iterator([0],[v],report_edges=True, labels=True))
            lcosets.append(braid_from_path(path, gens))
    return lcosets

def build_stab_gens(Gr, gens, pred=braid_in_img):
    lcosets = build_lcosets(Gr,gens)
    return gens_from_lcosets(gens,lcosets,pred=braid_in_img)

# Some basic group theory helpers

def conj(g,h):
    return h * g * h**(-1)

def commutator(g,h):
    return g*h*(g**(-1)) * (h**(-1))

def all_commutators(lst):
    return [commutator(g,h) for (g,h) in it.combinations(lst,2)]


# In[11]:


## Orbit Order with Braid Action ##

def find_orb_order(b,rho = None,N=1000):
    for i in range(1,N):
        if braid_in_img(b**i,rho):
            return i
    return None


# In[12]:


var('x')
matrix_from_scc(matrix([[1],[x]]))


# ## Hyperbolic Geometry

# In[13]:


def farey_next(t1,t2):
    (a,b,alpha) = t1
    (c,d,beta) = t2
    (x,y) = (a+c,b+d)
    k = gcd(x,y)
    return (x/k,y/k,N((alpha+beta)/2))

def farey_seqs(n):
    if n == 0:
        return [[(0,1,0),(1,0,pi)]]
    result = farey_seqs(n-1)
    farey_last = result[-1]
    next_farey = []
    farey_len = len(farey_last)

    for i in range(0,farey_len):
        next_farey.append(farey_last[i])
        if i < (farey_len-1):
            next_farey.append(farey_next(farey_last[i],farey_last[i+1]))

    return result + [next_farey]

def farey_angles(n):
    result = {}
    for (a,b,theta) in farey_seqs(n):
        result[(a,b)] = theta
    return result

def farey_seq(n):
    return [(a,b) for (a,b,_) in farey_seqs(n)[n-1]]

def pair_to_PD_point(t,H):
    UHP = H.UHP()
    PD = H.PD()

    (a,b) = t
    if b == 0:
        return PD.get_point(0.9999*I)
    else:
        return UHP.get_point(a/b+0.00001).to_model("PD")

def angle_to_PD_point(theta,H):
    PD = H.PD()
    return PD.get_point(0.999*exp(I*theta))

def farey_complex(H,n,mark_pts=[],make_dict=False):
    PD = H.PD()
    farey_pts_half = [[(a,b,angle_to_PD_point(theta,H)) for (a,b,theta) in lst] for lst in farey_seqs(n)]
    other_half = [[(-a,b,PD.get_point(pt.coordinates().conjugate())) for (a,b,pt) in lst] for lst in farey_pts_half]
    farey_pts = [lst1 + lst2 for (lst1,lst2) in zip(farey_pts_half,other_half)]

    mark_hyp_pts = [pair_to_PD_point(t,H) for t in mark_pts]
    farey_arcs = [H.PD().get_geodesic(p[2],q[2]) for lst in farey_pts for (p,q) in itertools.pairwise(lst)]
    #farey_arcs.append(H.PD().get_geodesic(farey_pts[-1],farey_pts[0]))
    farey_pts_dict = None
    if make_dict:
        farey_pts_dict = {}
        for (a,b,pt) in farey_pts[-1]:
            #z = pt.coordinates()
            #conj_pt = PD.get_point(z.conjugate())
            farey_pts_dict[(a,b)] = pt
            farey_pts_dict[(-a,-b)] = pt
            #farey_pts_dict[(-a,b)] = conj_pt
            #farey_pts_dict[(a,-b)] = conj_pt
        return (farey_arcs,mark_hyp_pts,farey_pts_dict)
    return (farey_arcs, mark_hyp_pts)

def display_farey(arcs,points,model="PD",arc_color="black"):
    new_arcs = [arc.to_model(model).plot(color=arc_color) for arc in arcs]
    new_pts = [pt.to_model(model).show(**opts) for (pt,opts) in points]
    res = new_arcs[0]
    for x in (new_arcs + new_pts)[1:]:
        res += x
    return res

def farey_with_curves(n=4,marked_curves=[],model="PD",preloaded=None):
    if preloaded:
        (farey_arcs, _, farey_pts_dict) = preloaded
    else:
        H = HyperbolicPlane()
        (farey_arcs, _, farey_pts_dict) = farey_complex(H,n,make_dict=True)
    return display_farey(farey_arcs, [(farey_pts_dict[t],opts) for (t,opts) in marked_curves],model=model)


# In[14]:


def farey_rho(rho,opts=None,preloaded=None):
    SCCs = [(s[0][0],s[1][0]) for s in scc_rep(rho)]
    if not preloaded:
        n = 0
        while not (all([s in farey_seq(n) for s in SCCs])):
            n += 1
        H = HyperbolicPlane()
        preloaded = farey_complex(H,n,make_dict=True)
    m = len(SCCs)
    color_dict = rainbow(m)
    SCCs_opts = [(s,{"color": c, "size": (i+1)*70, "alpha": 0.5}) for (s,c,i) in zip(SCCs,color_dict,range(0,m))]
    return farey_with_curves(marked_curves=SCCs_opts,preloaded=preloaded)


# In[15]:


import random
def farey_dist_tri_lst(lst):
    if len(lst) == 1:
        if lst[0] == 1:
            return 1
        else:
            return 2

    if lst[0] == 1:
        return farey_dist_tri_lst([1+lst[1]] + lst[2:])
    else:
        return farey_dist_tri_lst(lst[1:]) + 1

def farey_dist_infty(q):
    if q == Infinity:
        return 0
    elif q == 0 or q == 1:
        return 1
    elif 0 < q <= 1:
        continued_frac = list(abs(q).continued_fraction())
        return farey_dist_tri_lst(continued_frac)
    else:
        return farey_dist_infty(frac(abs(q)))

def mat_to_infty(r):
    if r == Infinity:
        return identity_matrix(2)
    p = r.numerator()
    q = r.denominator()
    _, a, b = xgcd(p,q)
    return matrix([[a,-b],[-q,p]])

def mat_on_rat(A,r):
    return curve_to_rat(A*rat_to_curve(r))

def farey_dist(p,q):
    A = mat_to_infty(p)
    return farey_dist_infty(mat_on_rat(A,q))

def farey_prod_dist(ps,qs):
    return sum([farey_dist(p,q) for (p,q) in zip(ps,qs)])

def braid_stretch(b,ps,qs):
    dist_before = farey_prod_dist(ps,qs)
    ps_new = apply_braid_rat(b,ps)
    qs_new = apply_braid_rat(b,qs)
    dist_new = farey_prod_dist(ps_new,qs_new)
    if dist_before == 0:
        return 0
    return dist_new/dist_before

def rand_ps(n,N=100):
    def rand_rat():
        a = random.randint(-N,N)
        b = random.randint(-N,N)
        if b == 0:
            return Infinity
        else:
            c = Rational(a/b)
            return c
    return [rand_rat() for i in range(0,n)]

def rand_braid_stretch(b,N=100):
    ps = rand_ps(b.strands(),N)
    qs = rand_ps(b.strands(),N)
    stretch = braid_stretch(b,ps,qs)
    return (ps, qs, stretch)

def many_rand_stretch(b,M=100,N=100,onlymax=True,printer=True):
    stretches = []
    for i in range(0,M):
        (ps, qs, stretch) = rand_braid_stretch(b,N)
        if not(onlymax) and printer:
            print("Trial ", i)
            print("  - ps = ", ps)
            print("  - qs = ", qs)
            print("  - stretch = ", stretch)
            print("")
        stretches.append(numerical_approx(stretch,digits=3))
    max_stretch = max(stretches)
    if printer:
        print("Max Stretch: ", max_stretch)
    return max_stretch

def matrix_stretch(mat,p):
    q = mat_on_rat(mat,p)
    return farey_dist(p,q)

def many_rand_mat_stretch(mat,M=100,N=100):
    ps = rand_ps(M,N)
    return max([matrix_stretch(mat,p) for p in ps])


# ## Action on Tangent

# In[16]:


# f is given as a list of values on generators gamma_1,...,gamma_n
# g is a group element in Fn
def cross_compute(f,g,func):
    F = g.parent()
    gens = g.parent().gens()
    assert len(gens) == len(f)
    lst = g.Tietze()
    l = len(lst)

    if l == 0:
        return f[0] - f[0]
    elif l == 1:
        i = lst[0]
        h = gens[abs(i)-1]**sgn(i)
        dude = f[abs(lst[0])-1]
        if lst[0] < 0:
            dude = -func(h, dude)
        return dude
    else:
        i = lst[0]
        h = gens[abs(i)-1]**sgn(i)
        rest = h**(-1)*g
        return cross_compute(f, h,func) + func(h, cross_compute(f,rest,func))

def cross_compute_act(f,b,rho):
    res = f
    new_rho = rho
    func = conj_act(new_rho)
    F_gens = FreeGroup(b.strands()).gens()
    B_gens = b.parent().gens()
    lst = b.Tietze()[::-1]
    for i in lst:
        tw = B_gens[abs(i)-1]**sgn(i)
        old_res = copy(res)
        res[abs(i)-1] = cross_compute(old_res, F_gens[abs(i)-1]*tw,func)
        res[abs(i)] = cross_compute(old_res, F_gens[abs(i)]*tw,func)
        new_rho = twist_rho_mat(tw,new_rho)
        func = conj_act(new_rho)
    return res


# In[17]:


def conj_act(rho):
    return lambda g,M: free_to_matrix(g,rho)*M*free_to_matrix(g**(-1),rho)

v1 = matrix([[1,0],[0,-1]])
v2 = matrix([[0,1],[0,0]])
v3 = matrix([[0,0],[1,0]])
z = matrix([[0,0],[0,0]])
n = 2
G = BraidGroup(2)
F = FreeGroup(2)
sigma = G.gens()[0]
gamma_1 = F.gens()[0]
gamma_2 = F.gens()[1]

V = QQbar^3
def sl2_mat_to_V(m):
    return (m[0,0],m[0,1],m[1,0])
def V_to_sl2_mat(v):
    e1 = matrix([[1,0],[0,-1]])
    e2 = matrix([[0,1],[0,0]])
    e3 = matrix([[0,0],[1,0]])
    return v[0]*e1+v[1]*e2+v[2]*e3

def do_it(pc,gam):
    func = conj_act(rho_std(n))
    comp = cross_compute(pc,gam*sigma**3,func)

    return sl2_mat_to_V(comp)

def bdry_deriv(v,rho):
    F = FreeGroup(len(rho))
    f = [V_to_sl2_mat(v[i:i+3]) for i in range(0,3*len(rho),3)]
    return cross_compute(f, prod(F.gens()), conj_act(rho))

def bdry_mat(rho):
    n = len(rho)
    V = QQ^(3*n)
    W = QQ^3
    lin =  linear_transformation(V,W,lambda v: sl2_mat_to_V(bdry_deriv(v,rho)))
    return lin.matrix().transpose()

def build_tangent_mat(g,rho=None,conj=matrix([[1,0],[0,1]])):
    n = g.strands()
    if rho == None:
        rho = rho_std(n)
    assert twist_rho_mat(g, rho=rho) == [conj*m*conj**(-1) for m in rho],"Not a Fixed Point"
    F = FreeGroup(n)
    gens = F.gens()
    z = matrix(2)

    func = conj_act(rho) 
    V = QQ^3
    def _mk_lst(v,j):
        lst = [z]*n
        lst[j] = V_to_sl2_mat(v)
        return lst
    comp_dict = dict([
        ((v,j),[conj**(-1)*m*conj for m in cross_compute_act(_mk_lst(v,j),g,rho)]) 
        for v in V.basis()
        for j in range(0,n)
    ])

    def _build_block(i,j):
        def _do_it_sub(v):
            lst = [z]*n
            lst[j] = V_to_sl2_mat(v)
            #comp = fox_it(gens[i]*g,lst,rho=rho)
            #comp = conj**(-1)*cross_compute(lst,gens[i]*g,func)*conj
            comp = comp_dict[v,j][i]
            return sl2_mat_to_V(comp)
        return linear_transformation(V,V,_do_it_sub).matrix().transpose()
    return block_matrix([[_build_block(i,j) for j in range(0,n)] for i in range(0,n)],subdivide=False)


def build_principal(v,rho):
    mat = V_to_sl2_mat(v)
    f = [mat-x*mat*x**(-1) for x in rho]
    return f

def build_principal_vec(v,rho):
    f = build_principal(v,rho)
    for s in [sl2_mat_to_V(m) for m in f]:
        yield from s

def trace_deriv(v,mat):
    var('t')
    return derivative(((t*V_to_sl2_mat(v)).exp() * mat).trace(),t).subs({t: 0})

def total_trace_deriv(n=None, rho=None):
    assert rho != None or n != None, "one must be defined"
    if rho == None:
        rho = rho_std(n)
    else:
        n = len(rho)
    V = QQ^3
    def block_guy(i):
        nonzero = [trace_deriv(x,rho[i]) for x in [(1,0,0),(0,1,0),(0,0,1)]]
        res = [0]*(3*i) + nonzero + [0]*(3*n-3*i-3)
        return res
    return matrix([block_guy(i) for i in range(0,n)])

def build_tangent_mod_princ(g,rho=None,conj=matrix([[1,0],[0,1]])):
    n = g.strands()
    if rho == None:
        rho = rho_std(n)
    mat = build_tangent_mat(g,rho=rho,conj=conj)
    V = QQ^(3*n)

    quots = [tuple(build_principal_vec(x,rho)) for x in [(1,0,0),(0,1,0),(0,0,1)]]
    W = V/V.subspace(quots)
    lft = W.lift_map()
    quot_map = W.quotient_map()
    return linear_transformation(W,W,lambda v: quot_map(mat*lft(v))).matrix().transpose()

def build_tangent_in_tr(g,rho=None,conj=matrix([[1,0],[0,1]])):
    n = g.strands()
    if rho == None:
        rho = rho_std(n)
    mat = build_tangent_mat(g,rho=rho,conj=conj)
    V = QQ^(3*n)
    W = total_trace_deriv(rho=rho).right_kernel()
    return linear_transformation(W,W,lambda v: mat*v).matrix().transpose()

def build_tangent_tr_princ(g,rho=None,conj=matrix([[1,0],[0,1]])):
    n = g.strands()
    if rho == None:
        rho = rho_std(n)
    mat = build_tangent_mat(g,rho=rho,conj=conj)
    V = QQ^(3*n)
    W = total_trace_deriv(rho=rho).right_kernel()
    #print(W)
    quots = [tuple(build_principal_vec(x,rho)) for x in [(1,0,0),(0,1,0),(0,0,1)]]
    Q = W/W.subspace(quots)
    #print(W.subspace(quots))
    lft = Q.lift_map()
    quot = Q.quotient_map()
    return linear_transformation(Q,Q,lambda v: quot(mat*lft(v))).matrix().transpose()

def build_tan_tr_tot_princ(g,rho=None,conj=matrix([[1,0],[0,1]])):
    n = g.strands()
    if rho == None:
        rho = rho_std(n)
    mat = build_tangent_mat(g,rho=rho,conj=conj)
    V = QQ^(3*n)
    trace_mat = total_trace_deriv(rho=rho)
    b_mat = bdry_mat(rho)
    lin = linear_transformation(V,QQ^(3+n),trace_mat.stack(b_mat),side='right')
    W = lin.kernel()
    quots = [tuple(build_principal_vec(x,rho)) for x in [(1,0,0),(0,1,0),(0,0,1)]]
    quots_space = V.subspace(quots)
    Q = W/W.intersection(quots_space)
    lft = Q.lift_map()
    quot = Q.quotient_map()
    return linear_transformation(Q,Q,lambda v: quot(mat*lft(v))).matrix().transpose()

def build_tan_tr_tot(g,rho=None):
    n = g.strands()
    if rho == None:
        rho = rho_std(n)
    mat = build_tangent_mat(g,rho=rho)
    V = QQ^(3*n)
    trace_mat = total_trace_deriv(rho=rho)
    b_mat = bdry_mat(rho)
    W = trace_mat.stack(b_mat).right_kernel()
    return linear_transformation(W,W,lambda v: mat*v).matrix().transpose()
