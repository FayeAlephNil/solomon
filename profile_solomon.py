"""
Flame-graph profiling for the solomon package.

Produces:
  profile_output/<scenario>.prof   — raw cProfile data
  profile_output/<scenario>.svg    — flame graph (open in a browser)

Usage:
  uv run python profile_solomon.py
  uv run python profile_solomon.py alternating_3 mod_rep_4   # subset
"""

import cProfile
import pstats
import io
import subprocess
import sys
from pathlib import Path

import solomon.examples as ex
import solomon.modular as smod
import solomon.utils as utils

OUT = Path("profile_output")
OUT.mkdir(exist_ok=True)


def run(label: str, fn):
    prof_path = OUT / f"{label}.prof"
    svg_path  = OUT / f"{label}.svg"

    print(f"\n{'='*60}")
    print(f"  Profiling: {label}")
    print(f"{'='*60}")

    pr = cProfile.Profile()
    pr.enable()
    fn()
    pr.disable()

    pr.dump_stats(str(prof_path))

    # Top-20 by cumulative time
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
    ps.print_stats(20)
    print(s.getvalue())

    # Flame graph via flameprof
    result = subprocess.run(
        ["uv", "tool", "run", "flameprof", str(prof_path), "-o", str(svg_path)],
        capture_output=True, text=True,
    )
    if result.returncode == 0:
        print(f"  Flame graph -> {svg_path}")
    else:
        print(f"  flameprof failed: {result.stderr.strip()}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Scenario helpers
# ---------------------------------------------------------------------------

def scenario_alternating_3():
    (_, _, orb) = ex.alternating_genus_1(3)
    orb.advance_until()

def scenario_alternating_4():
    (_, _, orb) = ex.alternating_genus_1(4)
    orb.advance_until()

def scenario_alternating_5():
    (_, _, orb) = ex.alternating_genus_1(5)
    orb.advance_until()

def scenario_hyper_ell_g2_k3():
    (_, _, orb) = ex.hyper_ell_fibration(2, 3)
    orb.advance_until()

def scenario_hyper_ell_g2_k4():
    (_, _, orb) = ex.hyper_ell_fibration(2, 4)
    orb.advance_until()

def scenario_mod_rep_3():
    (_, rep, orb) = ex.alternating_genus_1(3)
    prime_lst = utils.primes(50)
    mod_reps_orbs = smod.mod_rep_orbs(rep, orb, lst=prime_lst)
    mod_orbs = [o for (_, o) in mod_reps_orbs]
    smod.advance_all_orbs(mod_orbs)

def scenario_mod_rep_4():
    (_, rep, orb) = ex.alternating_genus_1(4)
    prime_lst = utils.primes(50)
    mod_reps_orbs = smod.mod_rep_orbs(rep, orb, lst=prime_lst)
    mod_orbs = [o for (_, o) in mod_reps_orbs]
    smod.advance_all_orbs(mod_orbs)


SCENARIOS = {
    "alternating_3":    scenario_alternating_3,
    "alternating_4":    scenario_alternating_4,
    "alternating_5":    scenario_alternating_5,
    "hyper_ell_g2_k3":  scenario_hyper_ell_g2_k3,
    "hyper_ell_g2_k4":  scenario_hyper_ell_g2_k4,
    "mod_rep_3":        scenario_mod_rep_3,
    "mod_rep_4":        scenario_mod_rep_4,
}

if __name__ == "__main__":
    if len(sys.argv) > 1:
        requested = sys.argv[1:]
        unknown = [r for r in requested if r not in SCENARIOS]
        if unknown:
            print(f"Unknown scenario(s): {unknown}")
            print(f"Available: {list(SCENARIOS)}")
            sys.exit(1)
        to_run = [(l, SCENARIOS[l]) for l in requested]
    else:
        to_run = list(SCENARIOS.items())

    for label, fn in to_run:
        run(label, fn)

    print(f"\nAll done. Open SVGs in {OUT}/ to view flame graphs.")
