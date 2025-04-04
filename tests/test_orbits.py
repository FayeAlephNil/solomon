import solomon
import solomon.examples as ex
import solomon.modular as smod
import solomon.utils as utils

### Integration Tests ###

def test_elliptic_disk_rep_3():
    (S, rep, orb) = ex.alternating_genus_1(3)
    orb.advance_until()
    assert len(orb.digraph.nodes) == 8
    assert len(orb.digraph.edges) == 16

def test_elliptic_disk_rep_4():
    (S, rep, orb) = ex.alternating_genus_1(4)
    orb.advance_until()
    assert len(orb.digraph.nodes) == 27

def test_ell_disk_mod_rep_3():
    (S,rep,orb) = ex.alternating_genus_1(3)
    prime_lst  = utils.primes(50)
    mod_reps_orbs = smod.mod_rep_orbs(rep,orb, lst=prime_lst)
    mod_orbs = [orb for (_,orb) in mod_reps_orbs]
    to_check = [len(orb.digraph.nodes) == 8 for orb in smod.advance_all_orbs(mod_orbs)]
    assert all(to_check)

def test_ell_disk_mod_rep_4():
    (S,rep,orb) = ex.alternating_genus_1(4)
    prime_lst  = utils.primes(50)
    mod_reps_orbs = smod.mod_rep_orbs(rep,orb, lst=prime_lst)
    mod_orbs = [orb for (_,orb) in mod_reps_orbs]
    to_check = [len(orb.digraph.nodes) == 27 for orb in smod.advance_all_orbs(mod_orbs)]
    assert all(to_check)
