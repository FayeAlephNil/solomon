# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Solomon is a Python research project for studying **Lefschetz fibrations** via their monodromy representations. The core question is how the automorphism group of smooth Lefschetz fibrations maps to the smooth mapping class group of the base surface. Paper-specific computations live in `PaperComputations/`; the main interactive workspace is `LefschetzDisk.ipynb`.

## Commands

```bash
# Install in editable mode
pip install -e .

# Run all tests
pytest

# Run a single test
pytest tests/test_orbits.py::test_elliptic_disk_rep_3

# CLI for orbit graph generation
python cli.py -n 3 -p 2 --output out.graph
python cli.py -N 5 -P 10 --dir ./output/  # batch over punctures and primes
```

## Architecture

### Core Package (`src/solomon/`)

The package builds up a pipeline: surface groups → monodromy representations → orbit graphs.

1. **`surface_group.py`** — Defines `FreeGrp`, `FreeGrpElement` (Tietze word representation), `FreeGrpHom`, and `SurfaceGroup` (genus g, n punctures). The fundamental building block.

2. **`monodromy_rep.py`** — `MonodromyRep` wraps a surface automorphism (from `curver`) with support for reducing it to the matrix action on first homology (`homology_matrices()`), as well as support for mod-N reductions (`to_local_system_mod()`).

3. **`orbit_graph.py`** — `OrbitGraph` builds directed multigraphs (via NetworkX) where nodes are representation elements and edges are generator actions. Supports multiprocessing for large graphs. Note that since this is an orbit graph, if the group-action is well-defined then the number of in edges and the number of out-edges from a node should both be equal to the number of generators.

4. **`modular.py`** — Batch helpers (`mod_rep_orbs()`, `advance_all_orbs()`) that run orbit computations across a range of primes.

5. **`examples.py`** — Factory functions for named Lefschetz fibration examples (`alternating_genus_1()`, `basic_chain()`, etc.), all using `curver` for surface automorphisms.

6. **`utils.py`** — Conversions between `python-flint` and NumPy matrices, prime sieve, cofactor computation.

### Supporting Systems

- **GAP scripts** (`.gap` files in root and `PaperComputations/`): Group computations run separately in GAP; results are used alongside the Python code.
- **SageMath** (`PaperComputations/*.sage`): Additional computations for papers, run in a Sage environment.
- **`cli.py`**: Click-based CLI wrapping `modular.py` batch functions.

### Key Dependencies

- `curver` / `flipper`: Core surface automorphism and mapping class group computations
- `python-flint`: Arbitrary-precision arithmetic for matrix operations
- `networkx`: Orbit graph structure
- `gappy-system`: Python↔GAP bridge (installed from git)
