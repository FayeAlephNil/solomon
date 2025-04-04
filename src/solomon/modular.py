from .monodromy_rep import MonodromyRep
from .orbit_graph import OrbitGraph

def mod_rep_orb(rep,orb,N):
    gens = orb.gens
    rep_N = rep.to_local_system_mod(N)
    orb_N = OrbitGraph(gens, rep_N, MonodromyRep.action)
    return (rep_N,orb_N)

def mod_rep_orbs(rep,orb,lst=range(2,10)):
    return [mod_rep_orb(rep,orb,N) for N in lst]

def advance_all_orbs(orbs,printit=False):
    if not(printit):
        for orb in orbs:
            orb.advance_until()
    else:
        for p,orb in orbs:
            print("# Working on #" + str(p))
            orb.advance_until()
            print(str(p) + " gives " + str(len(orb.digraph.nodes)) + " nodes")
    return orbs
