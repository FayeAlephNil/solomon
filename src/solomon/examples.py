# TODO: Refactor these as Lefschetz Fibrations, when that class is ready

import curver
from .surface_group import FreeGrp
from .monodromy_rep import MonodromyRep
from .orbit_graph import OrbitGraph
from .surface_group import SurfaceGroup

def alternating_genus_1(m):
    S = curver.load(1,1)
    g1 = S('a_0')
    g2 = S('b_0')
    pi1 = SurfaceGroup(0,m+1)
    rep = MonodromyRep(pi1, ([g1,g2]*m)[0:m])
    orb = OrbitGraph(FreeGrp(m-1).gens, rep, MonodromyRep.action)

    return (S, rep,orb)

def basic_chain(g,n):
    S = curver.load(g,n)
    chain_bs = [S(f"b_{i}") for i in range(0,g)]
    chain_cs = [S(f"c_{i}") for i in range(0,g-1)]

    chain = chain_bs + chain_cs
    chain[::2] = chain_bs
    chain[1::2] = chain_cs
    return (S,chain)

def chain_and_reverse(g,n):
    (S, chain) = basic_chain(g,n)
    chain.append(S(f'a_{g-1}'))
    chain.insert(0,S('a_0'))

    num_disk_punctures = 4*g
    pi1 = SurfaceGroup(0,num_disk_punctures+1)
    rep = MonodromyRep(pi1, 2*(chain + chain[::-1])[:-1])
    orb = OrbitGraph(FreeGrp(num_disk_punctures-1).gens, rep, MonodromyRep.action)
    return (S, rep,orb)

def partial_chain(g,n):
    (S, chain) = basic_chain(g,n)
    chain.insert(0,S('a_0'))

    repeat_num = 2*(2*g+1)
    num_chain = 2*g

    pi1 = SurfaceGroup(0,num_chain*repeat_num + 1)
    rep = MonodromyRep(pi1, chain * repeat_num)
    orb = OrbitGraph(FreeGrp(num_chain*repeat_num-1).gens, rep, MonodromyRep.action)
    return (S,rep,orb)

def full_chain(g,n):
    (S, chain) = basic_chain(g,n)
    chain.append(S(f'a_{g-1}'))
    chain.insert(0,S('a_0'))

    repeat_num = 2*g+2
    num_chain = 2*g+1

    pi1 = SurfaceGroup(0,num_chain*repeat_num + 1)
    rep = MonodromyRep(pi1, chain * repeat_num)
    orb = OrbitGraph(FreeGrp(num_chain*repeat_num-1).gens, rep, MonodromyRep.action)
    return (S,rep,orb)
