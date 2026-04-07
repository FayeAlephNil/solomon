"""
Tests for solomons_desk/ethan.py — Fibration class.

Covers two modes:
  - No-Sage  (always runs): GAP export, _schreier_generators, ImportError on Sage methods
  - Sage-available (skipped if Sage not installed): sage_braid_group,
    sage_stabilizer, sage_stabilizer_subgroup
"""

import pytest
from unittest.mock import patch

try:
    import sage.all  # noqa: F401
    HAS_SAGE = True
except ImportError:
    HAS_SAGE = False


def _make_fibration():
    """Genus-1, 4-twist alternating fibration (27-node orbit)."""
    from ethan import Fibration
    f = Fibration(genus=1, punctures=1)
    for c in ['a_0', 'b_0', 'a_0', 'b_0']:
        f.add_twist(c)
    return f


# ── No-Sage tests (always run) ────────────────────────────────────────────

def test_ethan_import():
    """Fibration is importable with no Sage present."""
    from ethan import Fibration  # noqa: F401


def test_gap_generators():
    """gap_generators() returns a valid GAP script string."""
    f = _make_fibration()
    f.build_orbit()
    script = f.gap_generators()
    assert isinstance(script, str)
    assert 'SurfaceBraidFpGroup' in script
    assert 'stab_gens' in script
    assert 'Subgroup' in script


def test_gap_generators_file(tmp_path):
    """gap_generators_file() writes a non-empty file."""
    f = _make_fibration()
    f.build_orbit()
    out = tmp_path / 'stab.g'
    f.gap_generators_file(str(out))
    assert out.exists()
    assert out.stat().st_size > 0
    assert 'SurfaceBraidFpGroup' in out.read_text()


def test_schreier_generators_type_and_content():
    """_schreier_generators() returns non-empty FreeGrpElement words."""
    from solomon.surface_group import FreeGrpElement
    f = _make_fibration()
    f.build_orbit()
    raw = f._schreier_generators()
    assert len(raw) > 0
    for w in raw:
        assert isinstance(w, FreeGrpElement)
        assert len(w.Tietze()) > 0  # no identity words


def test_schreier_consistent_with_gap_generators():
    """_schreier_generators() and gap_generators() must agree on count."""
    f = _make_fibration()
    f.build_orbit()
    raw = f._schreier_generators()
    script = f.gap_generators()
    # Count 's[' occurrences per-line is fragile; count stab_gens list entries
    lines = [l.strip() for l in script.splitlines() if l.strip().startswith('s[')]
    assert len(lines) == len(raw)


def test_sage_methods_raise_importerror_without_sage():
    """sage_braid_group/stabilizer/subgroup raise ImportError when Sage is blocked."""
    f = _make_fibration()
    f.build_orbit()
    # Block sage.all regardless of whether Sage is actually installed
    blocked = {'sage.all': None}
    with patch.dict('sys.modules', blocked):
        with pytest.raises(ImportError):
            f.sage_braid_group()
        with pytest.raises(ImportError):
            f.sage_stabilizer()
        with pytest.raises(ImportError):
            f.sage_stabilizer_subgroup()
        with pytest.raises(ImportError):
            f.sage_normal_core()
        with pytest.raises(ImportError):
            f.sage_core_quotient()


# ── Sage-available tests ──────────────────────────────────────────────────

@pytest.mark.skipif(not HAS_SAGE, reason="Sage not available in this environment")
def test_sage_braid_group_strand_count():
    """sage_braid_group() returns B_4 for a 4-twist fibration."""
    f = _make_fibration()
    B = f.sage_braid_group()
    assert B.strands() == 4


@pytest.mark.skipif(not HAS_SAGE, reason="Sage not available in this environment")
def test_sage_braid_group_before_orbit():
    """sage_braid_group() is callable before build_orbit()."""
    f = _make_fibration()
    # No build_orbit() call — braid group only needs twist count
    B = f.sage_braid_group()
    assert B.strands() == 4


@pytest.mark.skipif(not HAS_SAGE, reason="Sage not available in this environment")
def test_sage_stabilizer_tietze_roundtrip():
    """Sage Braid.Tietze() must exactly match the Solomon Tietze list."""
    f = _make_fibration()
    f.build_orbit()
    raw = f._schreier_generators()
    stab = f.sage_stabilizer()
    assert len(stab) == len(raw)
    for sage_b, sol_w in zip(stab, raw):
        assert list(sage_b.Tietze()) == sol_w.Tietze()


@pytest.mark.skipif(not HAS_SAGE, reason="Sage not available in this environment")
def test_sage_stabilizer_elements_in_correct_group():
    """Every stabilizer element must belong to sage_braid_group()."""
    f = _make_fibration()
    f.build_orbit()
    B = f.sage_braid_group()
    for b in f.sage_stabilizer():
        assert b.parent() == B


@pytest.mark.skipif(not HAS_SAGE, reason="Sage not available in this environment")
def test_sage_stabilizer_subgroup_has_generators():
    """sage_stabilizer_subgroup() returns a subgroup with generators."""
    f = _make_fibration()
    f.build_orbit()
    H = f.sage_stabilizer_subgroup()
    assert len(H.gens()) > 0


@pytest.mark.skipif(not HAS_SAGE, reason="Sage not available in this environment")
def test_sage_normal_core_index_divisible_by_orbit():
    """[B : Core_B(H)] is divisible by the orbit size."""
    from sage.all import libgap
    f = _make_fibration()
    f.build_orbit()
    gap_B = libgap(f.sage_braid_group())
    core = f.sage_normal_core()
    index = int(gap_B.Index(core))
    assert index % f.orbit_size() == 0


@pytest.mark.skipif(not HAS_SAGE, reason="Sage not available in this environment")
def test_sage_core_quotient_order():
    """B / Core_B(H) is a non-trivial finite group whose order is divisible by orbit_size."""
    f = _make_fibration()
    f.build_orbit()
    Q = f.sage_core_quotient()
    order = int(Q.Order())
    assert order > 1
    assert order % f.orbit_size() == 0
