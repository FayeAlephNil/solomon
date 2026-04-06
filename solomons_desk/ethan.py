"""
solomons_desk/ethan.py
High-level wrapper around solomon for Lefschetz fibration computations.

Usage:
    from ethan import Fibration

    fib = Fibration(genus=1, punctures=1)
    fib.add_twist('a_0')
    fib.add_twist('b_0')
    fib.add_twist('a_0')
    fib.build_orbit(mod=5)
    fib.show_orbit()
    fib.gap_generators_file('stab.g')
"""

from __future__ import annotations

import sys
import os
from collections import deque

import curver
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Make sure solomon is importable whether or not it is pip-installed
_repo_root = os.path.join(os.path.dirname(__file__), '..')
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from solomon.monodromy_rep import MonodromyRep
from solomon.orbit_graph import OrbitGraph
from solomon.surface_group import SurfaceGroup, FreeGrp


class Fibration:
    """
    A Lefschetz fibration over the disk, defined by a sequence of Dehn twists.

    Parameters
    ----------
    genus : int
        Genus of the fiber surface.
    punctures : int
        Number of punctures on the fiber surface (use 1 for a closed surface
        with one boundary component, which is the most common case).

    Examples
    --------
    Standard alternating genus-1 fibration with 3 twists::

        fib = Fibration(1, 1)
        fib.add_twist('a_0')
        fib.add_twist('b_0')
        fib.add_twist('a_0')
        fib.build_orbit(mod=5)
        fib.show_orbit()
    """

    def __init__(self, genus: int, punctures: int):
        self.genus = genus
        self.punctures = punctures
        self.surface = curver.load(genus, punctures)

        self._twists: list = []          # ordered list of curver MappingClass objects
        self._curve_registry: dict = {}  # name → MappingClass
        self._rep: MonodromyRep | None = None
        self._orbit: OrbitGraph | None = None
        self._braid_gens: list | None = None

    # ------------------------------------------------------------------
    # Defining curves and twists
    # ------------------------------------------------------------------

    def define_curve(self, name: str, curve) -> None:
        """
        Register a named mapping class for use with add_twist.

        Parameters
        ----------
        name : str
            The name you want to use later in add_twist.
        curve : str or curver MappingClass
            Either a curver curve name string (e.g. 'a_0') or a
            curver MappingClass object (including compositions).

        Examples
        --------
        Alias a standard curve::

            fib.define_curve('alpha', 'a_0')

        Define a custom mapping class as a composition::

            S = fib.surface
            fib.define_curve('conj', S('a_0') * S('b_0') * S('a_0') ** -1)
        """
        if isinstance(curve, str):
            self._curve_registry[name] = self.surface(curve)
        else:
            self._curve_registry[name] = curve
        self._invalidate()

    def add_twist(self, curve) -> None:
        """
        Append a Dehn twist generator to the fibration.

        Parameters
        ----------
        curve : str or curver MappingClass
            A curver curve name ('a_0', 'b_0', 'c_0', …), a name
            registered with define_curve, or a curver MappingClass object.

        Note
        ----
        You need at least 2 twists before calling build_orbit.
        """
        if isinstance(curve, str):
            if curve in self._curve_registry:
                mc = self._curve_registry[curve]
            else:
                mc = self.surface(curve)
        else:
            mc = curve
        self._twists.append(mc)
        self._invalidate()

    def clear_twists(self) -> None:
        """Remove all twists and reset the orbit."""
        self._twists.clear()
        self._invalidate()

    # ------------------------------------------------------------------
    # Building the orbit graph
    # ------------------------------------------------------------------

    def _invalidate(self) -> None:
        self._rep = None
        self._orbit = None

    def _build_rep(self) -> None:
        n = len(self._twists)
        if n < 2:
            raise ValueError(
                f"Need at least 2 twists to have a non-trivial orbit "
                f"(currently {n})."
            )
        domain = SurfaceGroup(0, n + 1)
        self._rep = MonodromyRep(domain, list(self._twists))
        self._braid_gens = FreeGrp(n - 1).gens

    def build_orbit(self, mod=None, max_nodes: int | None = None) -> OrbitGraph:
        """
        Compute the orbit of the monodromy representation under the braid group.

        Parameters
        ----------
        mod : None, -1, or positive int
            None  — use curver MappingClass objects directly (may be infinite)
            -1    — convert to integer homology matrices (no modular arithmetic)
            N > 0 — reduce mod N; produces a finite orbit suitable for GAP export
        max_nodes : int or None
            Stop as soon as the graph reaches this many nodes. Useful as a
            safety net when the orbit may be very large or infinite.
            None (default) — run until the orbit is fully explored or
            until you interrupt the kernel.

        Returns
        -------
        OrbitGraph
            The orbit graph, fully advanced or truncated at max_nodes.
        """
        if self._rep is None:
            self._build_rep()

        if mod is None:
            rep = self._rep
        elif mod == -1:
            rep = self._rep.to_local_system()
        else:
            rep = self._rep.to_local_system_mod(int(mod))

        self._orbit = OrbitGraph(self._braid_gens, rep)

        if max_nodes is None:
            self._orbit.advance_until()
        else:
            last_size = -1
            this_size = len(self._orbit.digraph.nodes())
            while last_size != this_size and this_size < max_nodes:
                last_size = this_size
                self._orbit.advance()
                this_size = len(self._orbit.digraph.nodes())
            if this_size >= max_nodes:
                print(f"Stopped at {this_size} nodes (max_nodes={max_nodes}). "
                      "The orbit may be larger; increase max_nodes to explore further.")

        return self._orbit

    def orbit_size(self) -> int:
        """Number of nodes in the current orbit graph."""
        self._require_orbit()
        return len(self._orbit.digraph.nodes())

    # ------------------------------------------------------------------
    # Visualisation
    # ------------------------------------------------------------------

    def _require_orbit(self) -> None:
        if self._orbit is None:
            raise RuntimeError("Call build_orbit() first.")

    # Layouts that accept a seed and those that don't
    _SEEDED_LAYOUTS = {
        'spring':      nx.spring_layout,
        'kamada_kawai': nx.kamada_kawai_layout,   # no seed, but harmless to pass
        'spectral':    nx.spectral_layout,
        'shell':       nx.shell_layout,
        'circular':    nx.circular_layout,
        'random':      nx.random_layout,
    }

    def _compute_layout(self, G, layout: str, seed: int):
        fn = self._SEEDED_LAYOUTS.get(layout)
        if fn is None:
            raise ValueError(
                f"Unknown layout '{layout}'. "
                f"Choose from: {list(self._SEEDED_LAYOUTS)}"
            )
        # Only spring and random accept a seed kwarg
        if layout in ('spring', 'random'):
            return fn(G, seed=seed)
        return fn(G)

    def _draw(self, ax, remove_loops: bool = True,
              layout: str = 'spring', seed: int = 42) -> None:
        self._require_orbit()
        G = nx.convert_node_labels_to_integers(self._orbit.digraph)
        if remove_loops:
            G.remove_edges_from(list(nx.selfloop_edges(G)))
        if len(G.nodes()) == 0:
            ax.text(0.5, 0.5, "Empty orbit graph", transform=ax.transAxes,
                    ha='center', va='center')
            return

        pos = self._compute_layout(G, layout, seed)
        all_edges = list(G.edges(keys=True))

        # Collect unique generator keys in encounter order
        unique_keys = list(dict.fromkeys(k for _, _, k in all_edges))
        palette = plt.cm.tab10.colors
        color_map = {k: palette[i % 10] for i, k in enumerate(unique_keys)}

        # Draw nodes and labels
        nx.draw_networkx_nodes(G, pos, ax=ax, node_size=350,
                               node_color='steelblue', alpha=0.9)
        nx.draw_networkx_labels(G, pos, ax=ax, font_color='white',
                                font_size=8, font_weight='bold')

        # Draw edges grouped by generator key, with distinct arc radii so
        # parallel edges between the same pair of nodes are visually separated
        n_keys = max(len(unique_keys), 1)
        for i, key in enumerate(unique_keys):
            key_edges = [(u, v, k) for u, v, k in all_edges if k == key]
            rad = 0.1 + 0.2 * i / n_keys
            nx.draw_networkx_edges(
                G, pos, ax=ax,
                edgelist=key_edges,
                edge_color=color_map[key],
                arrows=True,
                arrowsize=18,
                connectionstyle=f'arc3,rad={rad}',
                min_source_margin=12,
                min_target_margin=12,
            )

        # Legend
        handles = [
            Line2D([0], [0], color=color_map[k], lw=2,
                   label=f"σ_{i+1}  ({k})")
            for i, k in enumerate(unique_keys)
        ]
        ax.legend(handles=handles, loc='best', fontsize=8)
        ax.set_title(
            f"Orbit graph — {len(G.nodes())} nodes, "
            f"genus {self.genus}, {len(self._twists)} twists",
            fontsize=10,
        )
        ax.axis('off')

    def show_orbit(self, remove_loops: bool = True, figsize: tuple = (10, 8),
                   layout: str = 'spring', seed: int = 42) -> None:
        """
        Display the orbit graph inline (Jupyter / interactive).

        Parameters
        ----------
        remove_loops : bool
            Hide self-loop edges (default True).
        figsize : tuple
            Matplotlib figure size in inches.
        layout : str
            Node layout algorithm. Options: 'spring' (default), 'kamada_kawai',
            'spectral', 'shell', 'circular', 'random'.
        seed : int
            Random seed for layouts that use one ('spring', 'random').
            Adjust to get a cleaner-looking arrangement.
        """
        fig, ax = plt.subplots(figsize=figsize)
        self._draw(ax, remove_loops=remove_loops, layout=layout, seed=seed)
        plt.tight_layout()
        plt.show()

    def save_orbit_png(self, filename: str, remove_loops: bool = True,
                       dpi: int = 150, layout: str = 'spring',
                       seed: int = 42) -> None:
        """
        Save the orbit graph as a PNG image.

        Parameters
        ----------
        filename : str
            Output file path (e.g. 'orbit.png').
        remove_loops : bool
            Hide self-loop edges (default True).
        dpi : int
            Image resolution.
        layout : str
            Node layout algorithm (same options as show_orbit).
        seed : int
            Random seed for stochastic layouts.
        """
        fig, ax = plt.subplots(figsize=(10, 8))
        self._draw(ax, remove_loops=remove_loops, layout=layout, seed=seed)
        plt.tight_layout()
        fig.savefig(filename, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        print(f"Orbit graph saved to: {filename}")

    def save_orbit_dot(self, filename: str) -> None:
        """
        Save the orbit graph as a Graphviz DOT file.

        The file can be rendered externally with Graphviz::

            dot -Tpng orbit.dot -o orbit.png

        or loaded into other graph tools.
        """
        self._require_orbit()
        G = nx.convert_node_labels_to_integers(self._orbit.digraph)
        nx.drawing.nx_pydot.write_dot(G, filename)
        print(f"Orbit graph (DOT) saved to: {filename}")

    # ------------------------------------------------------------------
    # GAP export
    # ------------------------------------------------------------------

    def gap_generators(self) -> str:
        """
        Compute stabilizer generators of the orbit as braid group words,
        using the Schreier lemma (mirrors LefschetzDisk.ipynb build_stab_gens).

        The orbit must have been built with build_orbit(mod=N) first — the
        mod-N reduction is what makes the orbit finite; the returned generators
        are unreduced words in the braid group B_m, NOT matrices.

        Returns
        -------
        str
            A GAP script that defines a BraidGroup F, lists the stabilizer
            generators, and forms the subgroup H.
        """
        self._require_orbit()

        G = self._orbit.digraph
        start = self._orbit.start
        gen_map = {str(g): g for g in self._braid_gens}
        identity_word = self._braid_gens[0] ** 0

        # --- Step 1: BFS spanning tree → coset representatives ---
        # coset_rep[v] is a FreeGrpElement whose braid action sends `start` to v
        coset_rep = {start: identity_word}
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v, edge_dict in G[u].items():
                if v not in coset_rep:
                    key = next(iter(edge_dict))
                    coset_rep[v] = coset_rep[u] * gen_map[key]
                    queue.append(v)

        # --- Step 2: Schreier generators ---
        # For each edge u --g--> v, the element t_u * g * t_v^{-1} lies in
        # the stabilizer of `start`.
        stab_words = []
        for u, v, key in G.edges(keys=True):
            if u not in coset_rep or v not in coset_rep:
                continue  # shouldn't happen for a complete orbit
            g = gen_map[key]
            word = coset_rep[u] * g * (coset_rep[v] ** -1)
            if word.Tietze():  # skip identity
                stab_words.append(word)

        # --- Step 3: deduplicate by Tietze list ---
        seen: list = []
        for w in stab_words:
            if not any(w.Tietze() == s.Tietze() for s in seen):
                seen.append(w)

        # --- Step 4: Tietze → GAP braid word string ---
        def tietze_to_gap(tietze: list) -> str:
            if not tietze:
                return "Identity(F)"
            terms = []
            for i in tietze:
                if i > 0:
                    terms.append(f"s[{i}]")
                else:
                    terms.append(f"s[{-i}]^(-1)")
            return "*".join(terms)

        # --- Step 5: assemble GAP script ---
        m = len(self._braid_gens)
        mod_val = start.homology_mod_rep
        orbit_desc = f"mod {mod_val}" if mod_val != -1 else "full (curver)"

        lines = [tietze_to_gap(w.Tietze()) for w in seen]
        body = ",\n  ".join(lines) if lines else "# (trivial stabilizer)"

        script = (
            f"# Stabilizer of the monodromy representation in B_{m+1}\n"
            f"# Fibration: genus {self.genus}, {len(self._twists)} twists\n"
            f"# Orbit: {orbit_desc}  ({self.orbit_size()} nodes)\n"
            f"#\n"
            f"# s[i] denotes the i-th standard Artin generator sigma_i\n"
            f"LoadPackage(\"ACE\");\n"
            f"TCENUM := ACETCENUM;;\n"
            f"LoadPackage(\"FR\");\n"
            f"F := SurfaceBraidFpGroup({m+1},0,1);\n"
            f"s := GeneratorsOfGroup(F);\n"
            f"stab_gens := [\n"
            f"  {body}\n"
            f"];\n"
            f"H := Subgroup(F, stab_gens);\n"
        )
        return script

    def gap_generators_file(self, filename: str) -> None:
        """Write the GAP stabilizer script to a file."""
        code = self.gap_generators()
        with open(filename, 'w') as f:
            f.write(code)
        print(f"GAP stabilizer generators written to: {filename}")

    # ------------------------------------------------------------------
    # Misc
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        n = len(self._twists)
        return (f"Fibration(genus={self.genus}, punctures={self.punctures}, "
                f"twists={n})")
