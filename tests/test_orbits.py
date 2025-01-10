import solomon
import solomon.examples

### Integration Tests ###

def test_elliptic_disk_rep_3():
    (rep, orb) = solomon.examples.alternating_genus_1(3)
    orb.advance_until()
    assert len(orb.vertices()) == 8
    assert len(orb.edges()) == 16

def test_elliptic_disk_rep_4():
    (rep, orb) = solomon.examples.alternating_genus_1(4)
    orb.advance_until()
    assert len(orb.vertices()) == 27

