#!/usr/bin/env python
# coding: utf-8

# In[ ]:





# In[1]:


import itertools
n = 12
G = BraidGroup(n,'s')
G_gens = G.gens()
from operator import mul
# from functools import reduce # python3 compatibility
b = reduce(mul, G.gens(), G_gens[0]**0)
garside_to_prod = [b] + [reduce(mul,G_gens[:-i], G_gens[0]**0) for i in range(1,len(G_gens))]
garside = reduce(mul,garside_to_prod,G_gens[0]**0)
garside
meow = [reduce(mul,list(G_gens[:-i])[::-1], G_gens[0]**0) for i in range(1,len(G_gens))][::-1]+[reduce(mul, G.gens()[::-1], G_gens[0]**0)]
meow_braid=reduce(mul,meow,G_gens[0]**0)
r = reduce(mul,G.gens(),G_gens[0]**0)


# In[2]:


A = matrix([[1,-1],[0,1]])
B = matrix([[1,0],[1,1]])
rho_1 = [A,B,B,A,A,B,B,A]
G = BraidGroup(8,'s')
G_gens = G.gens()
rho_twist = twist_rho_mat(G_gens[1]*G_gens[2]*G_gens[1]^(-1) * G_gens[0]^2 * G_gens[5]*G_gens[6]*G_gens[5]^(-1),rho=rho_1)
scc_rep(rho_twist)


# In[ ]:


twist_rho_mat(r^(-1))


# In[42]:


twist_rho_mat(garside*r)


# In[55]:


from IPython.display import display
g = G_gens[0]
A = rho_std(2)[0]
B = rho_std(2)[1]
X = A
Y = matrix([[2,-1],[1,0]])
Z = B
conj = matrix([[1,-1],[1,0]])
rho_weird = [X,Y,Z]*4
(A*B*A)*(A*B*A)
(X*Y*Z)*(X*Y*Z)


# In[47]:


twist_rho(garside)


# In[5]:


BB = BraidGroup(2,'s')
F2 = FreeGroup(2)
F2_gens = F2.gens()
BB_gens = BB.gens()
scc_rep(twist_rho(BB_gens[0]**(-1)))
twist_rho(G_gens[0])


# In[63]:


many_sl2 = lots_of_sl2(10)


# In[31]:


idmat = matrix([[1,0],[0,1]])
alpha = matrix([[1],[0]])
C3 = matrix([[0,-1],[1,-1]])
C4 = matrix([[0,-1],[1,0]])
C6 = matrix([[1,-1],[1,0]])

matrix_from_scc(C3^2*alpha) * (matrix_from_scc(alpha)*matrix_from_scc(C3*alpha)*matrix_from_scc(C3^2*alpha))^(-4*d)


# In[35]:


w = matrix([[2],[1]])
matrix_from_scc(C6^2*w)*(matrix_from_scc(w)*matrix_from_scc(C6*w)*matrix_from_scc(C6^2*w))^(-4*d)


# In[19]:


var('p', 'q','r','s')
scc_1 = matrix([[p],[q]])
scc_2 = matrix([[r],[s]])
g = matrix([[0,-1],[1,-1]])
M = matrix_from_scc(scc_1) * matrix_from_scc(g * scc_1) * matrix_from_scc(g**2 * scc_1)
Npq = var('Npq')
M.trace()


# In[22]:


var('N')
expand(((N+2)*N -2))


# In[46]:


factor(p^4 + 2*p^2*q^2 + q^4-2)


# In[32]:


factor(p^4 - 2*p^3*q + 3*p^2*q^2 - 2*p*q^3 + q^4 - 2*p^2 + 2*p*q - 2*q^2)


# In[40]:


factor(p^4 + 2*p^3*q + 3*p^2*q^2 + 2*p*q^3 + q^4 + 2*p^2 + 2*p*q + 2*q^2)


# In[23]:


R.<xx,yy,zz,ww,N1N,N2N> = ZZ[]
S.<x,y,z,w,N1,N2> = R.quotient((xx**2 - xx*yy + yy**2 - N1N, zz**2 - zz*ww + ww**2 - N2N))
f = (a^2 - a*b + b^2 + 1)*(a^2 - a*b + b^2 - 2)^2*(a^2 - a*b + b^2 - 3)*Npr^2
floor(0.5)


# In[20]:


sccs = [scc_1,scc_2]
list_sccs = [g**(floor(i/2)) * sccs[i%2] for i in range(0,6)]
list_mats = [matrix_from_scc(scc) for scc in list_sccs]
our_mat = reduce(mul,list_mats,matrix([[1,0],[0,1]]))

#Mprime_lst = [matrix_from_scc(g**i * scc_1) * matrix_from_scc(g**i * scc_2) for i in range(0,3)]
#Mprime = reduce(mul,Mprime_lst, matrix([[1,0],[0,1]]))
#factor(Mprime.trace()**2 - 4)#.subs({p: xx, q: yy, r: zz, s: ww})
our_mat.trace().subs({s:0})


# In[27]:


@interact
def _(s=slider([-4 .. 4]), k=slider([0 .. 6]), r=slider([5..10])):
  var('x','y','z')

  F = our_mat.trace()(p=x,q=y,r=z,s=s)
  V = implicit_plot3d(F,(x,-r,r), (y,-r,r), (z,-r,r),
                     plot_points=120, color='red', smooth=True, countour=k)  
  show(V)


# In[32]:


x,y=var('x','y')
(our_mat.subs({r:1,s:-4}).trace()**2-4).factor().subs({p:x,q:y})
factor(p*q*r^2 - p^2*r*s - p*q*r*s + q^2*r*s + p^2*s^2 - p*q*s^2)


# In[138]:


(our_mat.trace()**2-4).factor()
((our_mat.subs({r:3,s:1}).trace())**2-4).factor()


# In[116]:


g2 = matrix([[0,-1],[1,0]])
list_sccs_2 = [g2**(floor(i/2)) * sccs[i%2] for i in range(0,4)]
list_mats_2 = [matrix_from_scc(scc) for scc in list_sccs_2]
our_mat_2 = reduce(mul,list_mats_2,matrix([[1,0],[0,1]]))


# In[119]:


g3 = matrix([[-1,-1],[1,0]])
list_sccs_3 = [g3**(floor(i/2)) * sccs[i%2] for i in range(0,6)]
list_mats_3 = [matrix_from_scc(scc) for scc in list_sccs_3]
our_mat_3 = reduce(mul,list_mats_3,matrix([[1,0],[0,1]]))


# In[15]:


alpha = matrix([[1],[0]])
beta = matrix([[0],[1]])
gamma = matrix([[1],[1]])
delta = matrix([[-1],[1]])
Ta = matrix_from_scc(alpha)
Tb = matrix_from_scc(beta)
Tc = matrix_from_scc(gamma)
Td = matrix_from_scc(delta)


# In[40]:


((Ta*Tc*Tb)*(Tb*Td*Ta))**2


# In[43]:


(Tc*Tb*Tb)


# In[47]:


Tc*Tb*Tb


# # Setup Cells

# ## Braid Group Action Setup

# In[2]:


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


# In[3]:


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


# In[4]:


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

# In[5]:


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


# In[6]:


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

# In[7]:


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



# In[8]:


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

# In[9]:


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


# # Example Graphs

# In[162]:


A = matrix([[1,-1],[0,1]])
B = matrix([[1,0],[1,1]])
hit_it = [A,B,B,A,A]
(Gr, gens4_weird, reps4_weird, timeout, comp_verts, _) = construct_orbit_hood(4,N=3000,M=3000,breadth=True, rho = [A,B,A,B])
Gr.remove_loops()
len(Gr.vertices())
display_orbit(Gr,gens4_weird)


# In[64]:


reps_enum = list(enumerate(reps4_weird))
S27 = SymmetricGroup(27)
S_gens = [S27([reps4_weird.index(twist_rho_mat(b,rho=p))+1 for (i,p) in reps_enum]) for b in gens4_weird]


# In[ ]:


ABAB
ABDA
ACBB
AACB


# In[63]:


S_gens[0] * S_gens_inv[0]


# In[43]:


scc_rep(reps4_weird[5])


# In[48]:


curves = set()
for r in reps4_weird:
    curves.update([(v[0,0],v[1,0]) for v in scc_rep(r)])
curves


# In[40]:


H = S27.subgroup(S_gens)
H.order()
H.blocks_all()


# In[31]:


reps4_weird.index(twist_rho_mat(gens4_weird[1],rho=reps4_weird[0]))


# In[17]:


C = matrix_from_scc(matrix([[1],[1]]))
D = matrix_from_scc(matrix([[-1],[1]]))
scc_rep(twist_rho_mat(BraidGroup(4).gens()[1], rho=[A,A,B,B]))


# In[22]:


C


# In[24]:


Gr = construct_orbit_hood(3,N=1000,M=1000,rho=[A,C,B])[0]
Gr.remove_loops()
len(Gr.vertices())


# In[64]:


import itertools
cands = [''.join(['a'] + list(x)) for x in itertools.product(['a','b'], repeat=4)]
cands


# In[96]:


gen_five = BraidGroup(5).gens()
scc_rep(twist_rho(gen_five[1]*gen_five[0],[A,B,B,B,B]))


# In[151]:


(Gr, gens, reps, timeout,comp_verts,_) = construct_orbit_hood(4,M=50,breadth=False,one_index_gen_names=True)
if timeout:
    print("Ran Out of Time")
#print(len(Gr.vertices()))
import sage.graphs.spanning_tree as spanning_tree
print(len(Gr.edges()))
print(len(spanning_tree.boruvka(Gr.to_undirected(),by_weight=False)))
Gr = Gr.subgraph([0,3,4,7,8,20,1,21,6])
print(len(Gr.edges()))
Gr.remove_loops()
nice_gr = display_orbit(Gr, gens,edge_labels=True,layout='spring',iterations=10000,heights={
    0: [1,20,0],
    1: [7,8,6],
    2: [3,21,4]
},save_pos=True)

#stab_gens = build_stab_gens(Gr,gens)
#identity = gens[0]*gens[0]**(-1)
#stab_gens = [s for s in stab_gens if s != identity]
#(huge_Gr, _, rreps, timeout, comp_verts2) = construct_huge_orbit(4,reps,gens=stab_gens)
#huge_Gr.remove_loops()
#display_orbit(huge_Gr, stab_gens, edge_labels=False).show(figsize=[2*6.4,2*4.8])


# In[153]:


nice_gr.show(figsize=[6.4,4.8])
Gr.get_pos()


# In[164]:


pos = Gr.get_pos()
pos[20] = [1.7,-0.2]
pos[8] = [2.08, 0.45]
#pos[7] = [1.48, 0.48]
#pos[6] = [2.68, 0.6]
pos[21] = [2.304, -0.65]
display_orbit(Gr, gens,edge_labels=True,pos=pos)


# ## Since we know that Core(G_4) = <<sigma^3>>
# stab_gens = [s for s in stab_gens if not (s.is_conjugated(gens[0]**3))]
# 
# def only_one(e,g = stab_gens[13]):
#     return e[2] == str(g)
# all_graphs = [huge_Gr.subgraph(edge_property = lambda e: only_one(e,g)) for g in stab_gens]
# #for gr in all_graphs:
#     #display_orbit(gr, stab_gens, edge_labels=False).show(figsize=[1*6.4,1*4.8])
# display_orbit(all_graphs[38], stab_gens, edge_labels=False).show(figsize=[2*6.4,2*4.8])

# In[64]:


(Gr,gens,reps,_,comp_verts,_) = construct_orbit_hood(5,M=400,N=400,breadth=False,rho=[A,A,A,A,B])
Gr.remove_loops()
#print(Gr.edges())
#display_orbit(Gr, gens).show()


# In[ ]:





# In[65]:


len(Gr.vertices())


# In[252]:


print([k for k,oth in enumerate(reps) if check_conj(reps[0],oth)[0]])
print([k for k,oth in enumerate(reps) if check_conj(reps[2],oth)[0]])
print([k for k,oth in enumerate(reps) if check_conj(reps[6],oth)[0]])
print([k for k,oth in enumerate(reps) if check_conj(reps[4],oth)[0]])
print([k for k,oth in enumerate(reps) if check_conj(reps[5],oth)[0]])
print([k for k,oth in enumerate(reps) if check_conj(reps[7],oth)[0]])


# In[278]:


_,conj = check_conj(reps[0],reps[3])
conj


# In[291]:


ell = [conj*m*conj**(-1) for m in twist_rho(BraidGroup(3).gens()[1], rho=reps[3])]
phi = twist_rho(BraidGroup(3).gens()[1],rho=reps[0])
print(ell)
scc_rep(twist_rho(BraidGroup(3).gens()[1], rho=reps[3]))


# In[300]:


check_conj(reps[2],reps[6])


# In[ ]:


cc = matrix([[0, a],[1, b]])
cc = 


# In[269]:


scc_rep(reps[7])


# In[261]:


check_conj(reps[4],reps[7])


# In[236]:


twist_rho_mat(BraidGroup(3).gens()[0]*BraidGroup(3).gens()[1])


# In[230]:


twist_rho(BraidGroup(3).gens()[1],rho=rho_std(3))


# In[232]:


reps[4]


# In[231]:


reps[4]


# In[229]:


var('a b')
conj = matrix([[0,1],[1,-2]])
conj*vector([0,1])


# In[26]:


(Gr, gens, reps4, timeout,comp_verts,_) = construct_orbit_hood(3,M=50,breadth=False)
if timeout:
    print("Ran Out of Time")
#print(len(Gr.vertices()))
#Gr = Gr.subgraph([0,3,4,7,8,20,1,21,6])
Gr.remove_loops()
#display_orbit(Gr, gens).show(figsize=[1.5*6.4,1.5*4.8])


# In[36]:


(Gr, gens, reps3, timeout,comp_verts,_) = construct_orbit_hood(4,N=2000,M=100,breadth=False)
if timeout:
    print("Ran Out of Time")
print(len(Gr.vertices()))
#Gr = Gr.subgraph([0,3,4,7,8,20,1,21,6])
#Gr.remove_loops()
#display_orbit(Gr, gens).show(figsize=[1.2*6.4,1.2*4.8])
#print(len(Gr.edges()))
(blah,blah2) = get_equi_antinodals(Gr, 0,gens)
#def long_prod(s,t,m):
#    res = s**0
#    for i in range(0,m):
#        if i % 2 == 0:
#            res = res*s
#        else:#
#            res = res*t
#    return res
#def find_m(s,t,N=200):
#    for m in range(1,N):
#        if long_prod(s,t,m) == long_prod(t,s,m):
#            return m
#    return None
#g = gens[0]**3
#h = gens[0]**(-1) * gens[1] * gens[0]
#k = gens[2] * gens[1] * gens[2]**(-1)
#print(find_m(g,k))
meow = [b.burau_matrix(reduced=True) for b in blah+blah2]
def find_weird(guys,N=10):
    e = guys[0] ** 0
    k = len(guys)
    guys_idxs = [i * j for j in range(1,k) for i in [-1,1]]
    combos = itertools.combinations_with_replacement(guys_idxs,N)
    for c in combos:
        g = from_Tietze(c,guys)
        if g != e:
            M = g.burau_matrix()
            I = identity_matrix(M.nrows())
            if M == I:
                return g
    return None
find_weird(blah+blah2,4)


# In[2]:


print(braid_in_img(blah[0]))
for x in blah + blah2:
    #print(path_from_braid(Gr, gens, 0, x))


# In[21]:


(Gr, gens5, reps5, timeout,comp_verts,_) = construct_orbit_hood(5,M=300,breadth=True)
if timeout:
    print("Ran Out of Time")
print(len(Gr.vertices()))
Gr.remove_loops()
#display_orbit(Gr, gens,completed_verts=comp_verts).show(figsize=[1*6.4,1*4.8])
#Gr.allow_loops(False)


# In[51]:





# In[26]:


rho_4_weird = nice_rep[0:4]
(Gr4_weird, gens4_weird, reps4_weird, timeout, comp_verts, _) = construct_orbit_hood(4,M=50,breadth=True, rho = rho_4_weird)
Gr4_weird.remove_loops()
display_orbit(Gr4_weird, gens4_weird,completed_verts=comp_verts).show(figsize=[5*6.4,5*4.8])


# In[33]:


meow = Gr4_weird.subgraph([1,5,13,37,9,19,21,3,8])
display_orbit(meow, gens4_weird).show(figsize=[5*6.4,5*4.8])


# In[34]:


(core_cand, not_core, not_core_ex) = find_core_tech(reps5,gens5,5)
core_cand


# In[9]:


#find_core_tech(reps5,gens5,6)


# # Scratchpad

# In[115]:


def central(gens):
    m = len(gens) + 1
    prod = gens[0]**0
    for g in gens:
        prod = prod * g
    return prod


# In[116]:


rep_free = [matrix([[1,-2],[0,1]]), matrix([[1,0],[2,1]])]
new_rep = rep_free


# ## Infinite Stuff

# In[17]:


a_path = Gr.all_paths(0,21,report_edges=True,labels=True)[0]
ell = braid_from_path(a_path,G_gens)


# In[416]:


Gr.longest_path(algorithm="heuristic").edges()


# In[453]:


c = G_gens[0]*G_gens[1]
braid_in_img(c**(-1) * G_gens[1]**3 * c)


# In[4]:


lcosets_3 = [G_gens[0]**0, G_gens[1],G_gens[0]*G_gens[1], G_gens[1]*G_gens[0]*G_gens[1],
             G_gens[0], G_gens[1]**2, G_gens[0]**2 * G_gens[1], G_gens[1]**2*G_gens[0]*G_gens[1]]
sub_gens = gens_from_lcosets(G_gens, lcosets_3)
#for g in sub_gens:
 #   g.plot().show()


# In[19]:


(Gr_cyc, gens_cyc, _,comp_verts) = construct_orbit_hood(5,[ell**(-1) *G_gens[3]**(-1)*G_gens[2]*ell], M=1000,breadth=False)
print(len(Gr_cyc.vertices()))
Gr_cyc.remove_loops()
#display_orbit(Gr_cyc, gens_cyc,edge_labels=False, completed_verts=comp_verts).show(figsize=[5*6.4,5*4.8])
Gr_cyc.allow_loops(False)
print("")


# In[19]:


from_ell = G_gens[3]**(-1) * G_gens[2]
special = ell**(-1) * from_ell * ell
#find_orb_order(special)


# In[11]:


ell


# In[18]:


twist_rho(ell)


# In[21]:


def some_niceties(m):
    return twist_rho(from_ell**m,twist_rho(ell))
def last_two(m):
    rho_lst = some_niceties(m)
    return rho_lst[-2]*rho_lst[-1]


# In[92]:


def check_special_sim(X,Y,N=1000):
    for i in range(-N,N):
        P = matrix([[1,0],[i,1]])
        if P*X*(P**(-1)) == Y:
            return (True,i)
    return False
check_special_sim(some_niceties(2)[2], some_niceties(12)[2])


# In[4]:


var('p')
var('q')
def make_ind_mat(x):
    return matrix([[x,-1],(x-1)**2,-(x-2)])
P = matrix([[p,-1],[(p-1)**2, -(p-2)]])
Q = matrix([[q,-1],[(q-1)**2, -(q-2)]])
simplify(expand(P**(-1)*Q*(P)))


# In[31]:


P_plus_one = matrix([[p+1,-1],[(p)**2, -(p-1)]])
P_minus_one = matrix([[p-1,-1],[(p-2)**2, -(p-1-2)]])
P_minus_two = matrix([[p-2,-1],[(p-3)**2, -(p-2-2)]])
factor(simplify(expand(P_plus_one**(-1)*P_minus_one * P_plus_one)))
factor(simplify(expand(P*P_plus_one * P**(-1))))


# In[6]:


expand((B)*P*(B**(-1)))


# ## Conjugacy Questions / Computations

# In[ ]:





# In[25]:


var('k')
K = matrix([[1,-1],[0,1]])
C = matrix([[0,-1],[1,2]])
D = matrix([[2,-1],[1,0]])
J = A*B*K
(A*B*K)*C*((A*B*K)**(-1))
(A*B*K)*A*(A*B*K)**(-1)
#conj_rho(J,twist_rho(G_gens[1],rho_std(3)))


# In[16]:


twist_rho(G_gens[0]**2)


# In[53]:


conj_lst = [matrix([[1,1],[-1,0]]),matrix([[0,-1],[1,1]])]
(helpers,alst) = test_conj_reps(reps, conj_lst,with_idx=True)


# In[54]:


len(alst)


# In[55]:


alst


# In[49]:


id(True)


# In[50]:


if []:



# In[66]:


reps[6]


# In[70]:


twist_rho(G_gens[2])


# In[76]:


reps[25]


# In[77]:


reps[19]


# In[95]:


for (j,r) in enumerate(reps):
    if r[1] == D and r[2] == D and r[3] == D:
        print(j)


# In[90]:


conj_rho(matrix([[0,-1],[1,-1]]), reps[19])


# In[84]:


D.is_similar(A,transformation=True)


# In[93]:


reps[13]


# In[14]:


meow = G_gens[1]*G_gens[0]*G_gens[1]**(-1)
meow2 = G_gens[1] ** (-1) *G_gens[0]*G_gens[1]
meow.plot().show()
meow2.plot().show()
twist_rho(meow)


# In[18]:


(G_gens[1]**3 * G_gens[0]**(-1) * G_gens[1] * G_gens[0]).plot()


# ## Specific Group Theory

# In[23]:


def construct_Gn(n):
    if n >= 5:
        return None
    (Gr, gens, _, _, _,BGroup) = construct_orbit_hood(n,M=100)
    #display_orbit(Gr,gens).show()
    stab_gens = build_stab_gens(Gr,gens)
    Gn = BGroup.subgroup(stab_gens)
    return (Gn,stab_gens, BGroup)



# In[13]:


(G3, G3_gens, B3) = construct_Gn(3)
B3_gens = B3.gens()
antinodals_G3 = [B3_gens[0]**3, 
                 B3_gens[1]**3, 
                 conj(B3_gens[0]**3, B3_gens[1]**(-1)),
                 conj(B3_gens[1]**3, B3_gens[0])]
equinodals_G3 = [conj(B3_gens[1], B3_gens[0]**(-1)),
                conj(B3_gens[0],B3_gens[1]),
                conj(B3_gens[1],B3_gens[1]**(-1) *B3_gens[0]),
                conj(B3_gens[0],B3_gens[0]*B3_gens[1]**(-1))
               ]
# Sanity check
assert(all([braid_in_img(b) for b in equinodals_G3]))
assert(all([braid_in_img(b) for b in antinodals_G3]))

equi_squares_G3 = [b**2 for b in equinodals_G3]
anti_squares_G3 = [b**2 for b in antinodals_G3]
G3
equi_squares_G3
#H3 = G3.quotient(equi_squares_G3 + anti_squares_G3)


# In[1]:


(G4, G4_gens, B4) = construct_Gn(4)
#G4_gens


# ## Chain Complexes

# In[14]:


from sage.topology.simplicial_complex import SimplicialComplex
from sage.modules.with_basis.morphism import ModuleMorphismFromFunction, ModuleMorphismByLinearity

(Gr3, gens3, _, _, _,BGroup) = construct_orbit_hood(3,M=70)

def build_chain_cx_edges(zero_chains, one_chains):
    zero_basis = zero_chains.basis()
    one_basis = one_chains.basis()
    def single_edge(edge):
        (s, t, name) = edge
        return zero_basis[t] - zero_basis[s]
    return ModuleMorphismByLinearity(domain=one_chains, codomain=zero_chains,on_basis = single_edge)

def create_relation_path(Gr, gens, v, i, j):
    if abs(i-j) > 1:
        lst = [(gens[i],1), (gens[j], 1), (gens[i],-1), (gens[j],-1)]
        return path_from_genlst_noswp(Gr, lst, v)
    else:
        lst1 = [(gens[i],1), (gens[j],1), (gens[i],1)]
        lst2 = [(gens[j],-1), (gens[i],-1), (gens[j],-1)]
        return path_from_genlst_noswp(Gr, lst1 + lst2, v)

def build_chain_cx_faces(Gr, gens, one_chains,two_chains):
    one_basis = one_chains.basis()
    two_basis = two_chains.basis()
    def chain_from_path(pth):
        tot = 0
        for (e, i) in pth:
            tot += i*one_basis[e]
        return tot

    def single_face(t):
        (v,i,j) = t
        pth = create_relation_path(Gr, gens, v, i, j)
        chain = chain_from_path(pth)
        return chain
    return ModuleMorphismByLinearity(domain=two_chains, codomain=one_chains,on_basis = single_face)


def build_chain_complex(n):
    (Gr, bgens, _, _, _,BGroup) = construct_orbit_hood(n,N=2000,M=100,breadth=True)
    zero_chains = FreeModule(ZZ, Gr.vertices())
    one_chains = FreeModule(ZZ, list(Gr.edges()))

    bgens_combs = itertools.combinations(range(0,len(bgens)),2)
    two_basis = [(v, i, j) for (i,j) in bgens_combs for v in Gr.vertices() ]
    two_chains = FreeModule(ZZ, two_basis)

    bdry_1_to_0 = build_chain_cx_edges(zero_chains,one_chains)
    bdry_2_to_1 = build_chain_cx_faces(Gr,bgens,one_chains,two_chains)

    return ChainComplex([two_chains, bdry_2_to_1.matrix(), one_chains, bdry_1_to_0.matrix(), zero_chains])

    if n >= 5:
        return None


# In[15]:


build_chain_complex(4).homology()


# In[46]:





# In[76]:





# In[ ]:


gens3[0].


# ## Finding Big Intersection

# In[96]:


(Gr, gens12, reps12, timeout,comp_verts,_) = construct_orbit_hood(12,N=2500,M=2500,breadth=True)


# In[97]:


def max_int(rho):
    return max(adj_int_numbers(rho))

lst_bigs = [(rho, max_int(rho)) for rho in reps12 if max_int(rho) > 2]
print(len(lst_bigs))
lst_bigs_better = [rho for rho in lst_bigs if A in rho]
#nice_rep = lst_bigs[0][0]
#meow = Gr.shortest_simple_paths(0,32,report_edges=True,labels=True)

#pth = next(meow)
#pth
#new_rep = nice_rep
#nice_rep
#pth
#reps5.index(nice_rep)
#pth
#lst_bigs


# In[98]:


def max_all_int(rho):
    return max(all_int_numbers(rho))
lst_bigs = [(rho, max_all_int(rho)) for rho in reps12 if max_all_int(rho) > 2]
print(len(lst_bigs))
nice_rep = lst_bigs[0][0]


# In[100]:


def max_eigen(m):
    return max([abs(l) for l in m.eigenvalues()])
def check_cand(rho):
    n = len(rho)
    def _check_cand(rho,j):
        record_it = max([max_eigen(prod(rho[:i])) for i in range(1,j)])
        first_bit = min([max_eigen(prod(rho[:i])) for i in range(1,j)])
        second_bit = min([max_eigen(prod(rho[j:i])) for i in range(j+1,n)])
        return (first_bit,second_bit,record_it)
    return [_check_cand(rho,j) for j in range(2,n-2)]
[check_cand(rho) for rho,_ in lst_bigs]


# In[185]:


nice_rho = twist_rho(gens13[3]*gens13[2]*gens13[1],lst_bigs[0][0])
scc_rep(nice_rho)
g = gens13[-2]
for i in range(0,3):
    print(max_int(twist_rho(g**i,nice_rho)))
    farey_rho(twist_rho(g**i,nice_rho),preloaded=preloaded).show()
#,preloaded=preloaded)


# ## Farey Fun

# In[57]:


H = HyperbolicPlane()
preloaded = farey_complex(H,6,make_dict=True)
lst_braids = [G_gens[3], G_gens[0],G_gens[0],G_gens[1]] + [G_gens[2]]*10
new_rho = rho_std(5)
many_rhos = [new_rho]
for g in lst_braids:
    many_rhos.append(twist_rho(g,rho=many_rhos[-1]))


# In[ ]:


farey_rho(new_rho,preloaded=preloaded).show()
for g in lst_braids:
    new_rho = twist_rho(g,rho=new_rho)
    farey_rho(new_rho,preloaded=preloaded).show()


# In[71]:


for r in many_rhos:
    for x in scc_rep(r):
        if x[1][0] != 0 and x[0][0] != 0:
            print(x[0][0]/x[1][0])
            print(farey_dist_infty(x[0][0]/x[1][0]))


# In[69]:


((-1/3).continued_fraction())[0]


# ## Tangent Computations

# In[96]:


G = BraidGroup(3)
G_gens = G.gens()
anti = G_gens[0]**3
boop = build_tangent_mat(G_gens[0]**3)
boop.eigenvalues()


# In[79]:


bop = build_tangent_mat(G_gens[0]**3)
bop.eigenvalues()


# In[78]:


boop.eigenvalues()


# In[113]:


m = 6*2
n = 2*m
G = BraidGroup(n)
G_gens = G.gens()
def twist_ij(i,j,gens):
    n = len(gens)+1
    i = i % n
    j = j % n
    if i > j:
        return twist_ij(j,i,gens)**(-1)
    if j == i + 1:
        return gens[i]
    conj = prod([gens[k] for k in range(i,j-1)])
    return conj*gens[j-1]*conj**(-1)
red_gens = [twist_ij(2*i,2*i+2,G_gens)**(-1) for i in range(0,m)]
blue_gens = [twist_ij(2*i+1,2*i+3,G_gens) for i in range(0,m)]
red_curve_twist = prod(red_gens[:-1])**m
blue_curve_twist = prod(blue_gens[:-1])**m
check_conj(twist_rho_mat(blue_curve_twist),rho_std(n))


# In[114]:


conj1 = matrix([[1,-m],[0,1]])
conj2 = matrix([[1,0],[-m,1]])
blue_mat = build_tan_tr_tot_princ(blue_curve_twist,conj=conj1**(-1))
red_mat = build_tan_tr_tot_princ(red_curve_twist,conj=conj2**(-1))


# In[105]:


conj1 = matrix([[1,-m],[0,1]])
conj2 = matrix([[1,0],[-m,1]])
get_ipython().run_line_magic('prun', 'blue_mat = build_tan_tr_tot_princ(blue_curve_twist,conj=conj1**(-1))')
red_mat = build_tan_tr_tot_princ(red_curve_twist,conj=conj2**(-1))


# In[39]:


conj1 = matrix([[1,-m],[0,1]])
conj2 = matrix([[1,0],[-m,1]])
get_ipython().run_line_magic('prun', 'blue_mat = build_tan_tr_tot_princ(blue_curve_twist,conj=conj1**(-1))')
red_mat = build_tan_tr_tot_princ(red_curve_twist,conj=conj2**(-1))


# In[ ]:





# In[75]:


B**6*A*B**(-6)


# # QuasiIsom Tests

# In[64]:


K = 8
collect_many = [many_rand_stretch(H_gens[0]**i, M=50, N=1000000, printer=False) for i in range(-K,K)]


# In[65]:


print(collect_many)
scatter_plot(many_logged)


# In[72]:


many_N = [many_rand_stretch(H_gens[0], M = 100, N = 2**n,printer=False) for n in range(0,25)]
scatter_plot(list(zip(range(0,20),many_N)))


# In[95]:


many_N


# 

# In[121]:


A = matrix([[1,-1],[0,1]])
B = matrix([[1,0],[1,1]])
fixed_qs = [4,0]
many_rand_ps = sum([rand_ps(5,N=5**i) for i in range(0,10)],[])
pss = [list(ps) for ps in itertools.combinations(many_rand_ps,2)]
calc_stretches = [numerical_approx(braid_stretch(H_gens[0], fixed_qs, ps)) for ps in pss]
mmax = max(calc_stretches)
print(mmax)
pss[calc_stretches.index(mmax)]


# In[122]:


apply_braid_rat(H_gens[0],[4,0])


# In[113]:


apply_braid_rat(H_gens[0],[Infinity,0])


# In[124]:


farey_dist(16/3,-1)


# In[15]:


apply_braid_rat(H_gens[0]**20,[-1,1])


# In[77]:


mmeow = []
keep = rand_ps(200,N=2**100000)
for p in keep:
    meow1 = [0,p]
    meow2 = [Infinity,p]
    mmeow.append(numerical_approx(braid_stretch(H_gens[0],meow1,meow2)))
meow_max = max(mmeow)
print(meow_max)
keep[mmeow.index(meow_max)]


# # Level 2 Graphs

# In[2]:


import networkx as nx

def clean_graph(file_name):
    nx_gr = nx.drawing.nx_pydot.read_dot(file_name)
    edge_list = list(nx_gr.edges(keys=True))
    edge_list_edit = [(u,v,k[0] + str(int(k[1])-1)) for (u,v,k) in edge_list]
    new_graph = DiGraph(edge_list_edit,loops=True,multiedges=True)
    N = len(new_graph.vertices())
    new_graph.relabel(dict([(str(i),i) for i in range(0,N+1)]))
    return new_graph

def build_gens(n,p,file_name,with_inv=True):
    if not(with_inv):
        print("Not implemented without inversion from OG graph")
        return
    gr = clean_graph(file_name)
    B = BraidGroup(n, 'a')
    B_gens = B.gens()
    lcos = build_lcosets(gr,[b**(-1) for b in B_gens])
    rcos = rcosets_from_lcosets(lcos)
    return gens_from_rcosets(B_gens, rcos, pred=lambda b: braid_in_img_p(b,p=p))


# In[13]:


B5 = BraidGroup(5,'a')
B5_gens = B5.gens()
the_graph = clean_graph('./output/gen_graphs_2s/5-2-out.graph')
the_graph.remove_loops()
nice_gr = display_orbit(the_graph, B5_gens)
nice_gr.show(figsize=[2*6.4,2*4.8]) 


# In[ ]:


ababaa_two_gens = set(build_gens(6,2,'./6_ababaa.dot'))


# In[ ]:


five_seven_gens = set(build_gens(5,7,'./five_seven.out'))


# In[ ]:


five_seven_gens


# In[ ]:


level_2_five_equianti = get_equi_antinodals(level_2_five_punct, 0, [b for b in B5_gens])


# In[45]:


#lst_cands = []
#many_braids = lots_of_braid(B5_gens,6)
#many_reps = [twist_rho(b) for b in many_braids]


# In[52]:


g = B5_gens[1] * B5_gens[0]**2 * B5_gens[3]
braid_in_img(g**(-1) * B5_gens[2]**1 * g)


# In[65]:


h = B5_gens[0]**2*B5_gens[2]*B5_gens[1]**2*B5_gens[2]*B5_gens[3]*B5_gens[2]*B5_gens[3]**(-1)*B5_gens[1]**(-2)*B5_gens[3]**(-1)*B5_gens[2]**(-1)*B5_gens[0]**(-2)


# In[67]:


braid_in_img_p(h)


# In[64]:


h


# In[6]:


level_2_five_punct.remove_loops()
level_2_five_punct.plot()


# In[28]:


A1 = matrix([[1,1],[0,1]])
A2 = matrix([[1,0],[1,1]])
A3 = matrix([[1,0],[0,1]])
A4 = matrix([[1,1],[1,0]])
A5 = matrix([[0,1],[1,1]])
A6 = matrix([[0,1],[1,0]])
mod_2_sl2 = [A1,A2,A3,A4,A5,A6]
pos_reps = [[a,b] for a in mod_2_sl2 for b in mod_2_sl2]
#[braid_in_img_p(BraidGroup(2).gens()[0]**3, rho=rep,p=2) for rep in pos_reps]


# In[14]:


#[braid_in_img((B5_gens[0]* B5_gens[1]*B5_gens[2])**4 ,rho=rep) for rep in reps5]


# # Finite Orders, weirdness, and graphs

# In[5]:


import sp4_torsion_work as finite_orders


# In[22]:


import math
math.prod([matrix_from_scc(matrix([[1],[-k]])) for k in range(0,9)])


# In[24]:


abs(1597-19075)


# In[56]:


rep = [A,A,A,A,B]
B5 = BraidGroup(5)
B5_gens  = B5.gens()
scc_rep(twist_rho_mat(B5_gens[3],rho=rep))
B


# In[68]:


for x in lots_of_braid(B5_gens,4):
    if any([p > 1 for p in all_int_numbers(twist_rho_mat(x,rho=rep))]):
        print(x)


# In[74]:


scc_rep(twist_rho_mat(B5_gens[1] * B5_gens[2] * B5_gens[3]**(-1),rho=rep))


# In[17]:


def lem_edge(n,startA=True,string=False):
    if n == 0:
        return [[]]
    A = matrix([[1,-1],[0,1]])
    B = matrix([[1,0],[1,1]])
    if string:
        A = 'A'
        B = 'B'
    lst = lem_edge(n-1,startA=False,string=string)
    if startA:
        return [[A] + x for x in lst]
    else:
        return [[A] + x for x in lst] + [[B] + x for x in lst]

lem_edge(4,string=True)


# In[44]:


def mat_to_hashable(m):
    return tuple(tuple(row) for row in m.rows())

def rep_to_hashable(rho):
    return tuple([mat_to_hashable(m) for m in rho])
n = 4
seen = set([])
(Gr_noconj,jens,reps_std,_,_,_) = construct_orbit_hood(n,M=400,N=400,breadth=False,rho=rho_std(n))
graphs=[]
new_stuff = []
for rep in lem_edge(n):
    if any([x != A for x in rep]):
        (Gr,gens,new_reps,_,_,_) = construct_orbit_hood_conj(n,M=400,N=400,breadth=False,rho=rep)
        L = len(Gr.vertices())
        if L < 200:
            graphs.append((Gr,gens))
            if all([not(check_conj(x,rep)[0]) for x in reps_std]):
                print(L)
                print(scc_rep(rep))
                new_stuff.append(rep)


# In[45]:


new_stuff


# In[37]:


display_orbit(*graphs[0])


# In[77]:


2**(3-1)


# In[29]:


B5 = BraidGroup(5)
B5_gens = B5.gens()
scc_rep(twist_rho_mat(B5_gens[1]*B5_gens[0],rho=[A,B,B,B,B]))


# # Finding Pants

# In[17]:





# In[18]:


scc_rep(rep)


# In[17]:


import itertools
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
        print(pants_idx)
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


# In[18]:


import itertools

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


# In[25]:


from operator import mul
import itertools

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
        print(d)


#build_special_rep(2)[12:]
bd_trace_pants(build_from_pants(example_pants(1).to_lst(),build_special_rep(1))["build"][0])


# In[26]:


2**2 + (-42)**2 + (-37)**2 - (-37)*(-42)*2


# In[67]:


def test_baby(j):
    lst = [build_special_rep(3)[i+4] for i in range(0,j+1)]
    return reduce(mul,lst,matrix([[1,0],[0,1]])).trace()
(build_special_rep(1)[4]*build_special_rep(1)[5]).trace()
test_baby(16)


# In[314]:


X = matrix([[7,-4],[9,-5]])
Y=X*B*A*A*A
print(Y)
Y*B*A*B*A*B*A


# # Finite Order Paper Check

# In[57]:


C4 = matrix([[0, -1], [1,0]])
C3 = matrix([[0,-1],[1,-1]])
C6 = matrix([[1,-1],[1,0]])
var('p q')
Tv = matrix_from_scc(matrix([[p],[q]]))
f4 = (Tv*C4*Tv*C4**(-1)).trace()
f3 = (Tv * C3 * Tv * C3**(-1) * C3**2 * Tv * C3**(-2)).trace()
f6 = (Tv * C6 * Tv * C6**(-1) * C6**2 * Tv * C6**(-2)).trace()


# In[76]:


N3 = p^2 - p*q + q^2
N4 = p^2 + q^2
f3_cand = -(N3 + 1)*(N3**2 + 2*N3 -2)
f4_cand = 2-N4**2
f6_cand = (N3 - 1)*(N3**2 -2 * N3 - 2)

assert (f3 - f3_cand).full_simplify() == 0
assert (f6 - f6_cand).full_simplify() == 0
assert (f4 - f4_cand).full_simplify() == 0


# In[26]:


from operator import mul

alpha = matrix([[1],[0]])
beta = matrix([[0],[1]])
gamma = matrix([[1],[1]])
delta = matrix([[1],[-1]])

A,B,C,D = tuple([matrix_from_scc(x) for x in [alpha,beta,gamma,delta]])

d = 1
n = 12*d
rho_n3 = [A,C,B,B,D,A]*(n//6)

G = BraidGroup(n)
r = reduce(mul,G.gens(),(G.gens())[0]**0)
check_conj(rho_n3,twist_rho_mat(r**3,rho_n3))


# In[ ]:





# In[ ]:





# # Legendre Family

# In[17]:


def check_conj_squares(rho1,rho2):
    vs = squares_lst(rho1)
    ws = squares_lst(rho2)
    if [v[0] for v in vs] != [w[0] for w in ws]:
        return False,None
    vecs1 = [v[1] for v in vs]
    vecs2 = [w[1] for w in ws]
    all_1 = all([int_number(v1,v2) == 1 for (v1, v2) in itertools.combinations(vecs1,2)])
    all_2 = all([int_number(v1,v2) == 1 for (v1, v2) in itertools.combinations(vecs2,2)])
    if all_1 != all_2:
        return False,None
    def is_counter(bs):
        matrix_from_scc(bs[1]) * bs[0] == bs[2]
    return is_counter(vecs1) == is_counter(vecs2),None

def squares_lst(rho):
    sgns = [m.eigenvalues()[0] for m in rho]
    def rationalize(m):
        a,b = m.right_eigenvectors()[0][1][0]
        d = lcm(a.denominator(), b.denominator())
        return matrix([[d*a],[d*b]])
    vecs = [rationalize(m) for m in rho]
    return list(zip(sgns,vecs))


# In[63]:


legendre_rho = [matrix([[1, 2], [0,1]]), matrix([[1, 0],[-2,1]]), matrix([[1,-2],[2,-3]])]
B3 = BraidGroup(3)


# In[ ]:





# In[64]:


twisted=twist_rho_mat(gens3[0]*gens3[1],rho=legendre_rho)
check_conj_squares(legendre_rho,twisted)


# In[65]:


(Gr, gens3, reps3_leg, _, comp_verts, _) = construct_orbit_hood_conj(3,N=30,M=30,breadth=True, rho = legendre_rho,conj_checker=check_conj_squares)
len(Gr.vertices())
display_orbit(Gr,gens3)


# In[68]:


matrix_from_scc(matrix([[1],[0]])).eigenvalues()


# In[33]:


hesse_rho= [matrix_from_scc(matrix([[2],[1]]))**3, 
            matrix_from_scc(matrix([[1],[1]]))**3, 
            matrix([[1,0],[1,1]])**3,
            matrix([[1,-1],[0,1]])**3]


# In[37]:


general_conj(build_farey_gens(3),hesse_rho)


# In[44]:


X*Z*Y


# In[35]:


[matrix([[1, 2], [0,1]]), matrix([[1, 0],[-2,1]]), matrix([[1,-2],[2,-3]])]


# In[18]:


def build_farey_gens(N,powwow=None):
    if powwow==None:
        powwow=N
    cusps = FareySymbol(Gamma(N)).cusps()
    def mat_from_cusp(c):
        S = matrix([[1,-powwow],[0,1]])
        if c == Infinity:
            return S
        else:
            s = c.numerator()
            t = c.denominator()
            g,p,q = xgcd(s,t)
            conj = matrix([[s,-q],[t,p]])
            return conj * S * conj**(-1)
    cusps.reverse()
    return [mat_from_cusp(c) for c in cusps]

build_farey_gens(2)


# In[261]:


(Gr, gens3, reps3_leg, _, comp_verts, _) = construct_orbit_hood_conj(6,N=200,M=10,breadth=True, rho = build_farey_gens(4),conj_checker=general_conj)


# In[260]:


comp_verts


# In[259]:


display_orbit(Gr.subgraph(Gr.neighbors(0)+[0]+Gr.neighbors(1)+Gr.neighbors(4)+Gr.neighbors(2)+Gr.neighbors(3)+Gr.neighbors(5)),gens3,edge_labels=True)


# In[144]:


reps3_leg[5]


# In[30]:


G = BraidGroup(4)
x = G.gens()[0]*G.gens()[1]*G.gens()[2]
y = G.gens()[1]*G.gens()[0]
(Gr_new_gens, gen_it, fun_reps, _, _, _) = construct_orbit_hood_conj(4,gens=[x,y],N=20,M=50,breadth=True, rho = build_farey_gens(3),conj_checker=general_conj)


# In[278]:


c = G.gens()[1] * G.gens()[2] * G.gens()[2]**(-1)*G.gens()[1]*G.gens()[2]
general_conj(fun_reps[10],twist_rho_mat(c**2,rho=fun_reps[10]))


# In[31]:


display_orbit(Gr_new_gens,gen_it)


# In[22]:


def general_scc(rho):
    def _do_one(M):
        a,b = M.right_eigenvectors()[0][1][0]
        d = lcm(a.denominator(), b.denominator())

        return matrix([[d*a],[d*b]])
    return [_do_one(M) for M in rho]

def general_conj(rho1,rho2,printer=False):
    return check_conj(rho1,rho2,make_scc=general_scc,printer=printer)


# In[76]:


scc_rep(rho_std(3))


# In[131]:


reduce(mul,build_farey_gens(5,powwow=5-1e-5),matrix([[1,0],[0,1]]))


# In[26]:


reduce(mul,[matrix([[-1,0],[0,-1]]) * x for x in build_farey_gens(5)],matrix([[1,0],[0,1]]))


# In[28]:


[x.trace() for x in build_farey_gens(4)]


# In[ ]:




