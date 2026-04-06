# Claude's Implementation Spec: Ethan Pesikoff Wrapper for Solomon

## Goal

Provide a clean, JupyterLab-friendly wrapper around solomon for Ethan Pesikoff (math PhD student, comfortable in JupyterLab/Mathematica, not a programmer). All files live in a self-contained subdirectory.

---

## Directory Name Suggestions (pick one)

1. **`fibration_lab`** — a workspace for fibration experiments
2. **`lefschetz_desk`** — Ethan's Lefschetz workstation
3. **`orbit_workshop`** — where orbit graphs are built
4. **`dehn_tools`** — Dehn twist computation tools

---

## Files to Create

```
<chosen-dir>/
├── ethan.py          # High-level wrapper module (~230 lines)
└── workspace.ipynb   # Starter notebook with worked examples
```

---

## `ethan.py` Design

### Key design principle: custom curves must be real mapping classes

Custom curves are defined as **compositions of curver `MappingClass` objects** (Dehn twists and their products), not homology matrices. This preserves the full mapping class structure needed by `MonodromyRep.twist()` and ensures correctness.

Ethan accesses the underlying curver surface via `fib.surface` and composes standard Dehn twists:
```python
fib = Fibration(1, 1)
S = fib.surface          # curver surface object
custom = S('a_0') * S('b_0') * S('a_0') ** -1   # conjugation
fib.define_curve('alpha_conj', custom)
fib.add_twist('alpha_conj')
```

---

### Class `Fibration`

```python
class Fibration:
    def __init__(self, genus, punctures)
    def define_curve(self, name, mapping_class)  # name a curver MappingClass
    def add_twist(self, curve)                   # str name or curver MappingClass
    def build_orbit(self, mod=None)              # → OrbitGraph, calls advance_until
    def orbit_size(self)                         # number of nodes in graph
    def show_orbit(self, figsize=(10,8), layout='spring')  # inline matplotlib
    def save_orbit_png(self, filename, dpi=150)  # PNG to disk
    def gap_generators(self, mod=None)           # → GAP-formatted string
    def gap_generators_file(self, filename, mod=None)      # write to .g file
```

#### `__init__`
```python
self.genus = genus
self.punctures = punctures
self.surface = curver.load(genus, punctures)   # curver surface
self._twists = []             # list of curver MappingClass objects
self._curve_registry = {}     # name → curver MappingClass
self._rep = None              # cached MonodromyRep
self._orbit = None            # cached OrbitGraph
self._braid_gens = None       # cached FreeGrp generators
```

#### `define_curve(name, mapping_class)`
Stores a named alias for a curver `MappingClass` (or composition thereof):
```python
self._curve_registry[name] = mapping_class
self._rep = self._orbit = None  # invalidate cache
```

#### `add_twist(curve)`
- `str`: look up in `_curve_registry` first, then `self.surface(curve)` for standard names (`'a_0'`, `'b_0'`, `'c_0'`, etc.)
- curver `MappingClass`: stored directly
- Invalidates `_rep` and `_orbit` caches

#### `_build_rep()` (internal)
Follows the pattern in `examples.py`:
```python
n = len(self._twists)
domain = SurfaceGroup(0, n + 1)                   # surface_group.py
self._rep = MonodromyRep(domain, self._twists)    # monodromy_rep.py:24
self._braid_gens = FreeGrp(n - 1).gens            # surface_group.py
```

#### `build_orbit(mod=None)`

Three modes controlled by `mod`:

| `mod` value | Behaviour |
|-------------|-----------|
| `None` (default) | Use curver `MappingClass` objects directly — no homology conversion |
| `-1` | Convert to integer homology matrices via `to_local_system()` — no modular arithmetic |
| positive int N | Convert to mod-N homology matrices via `to_local_system_mod(N)` — finite orbit |

```python
if self._rep is None:
    self._build_rep()
if mod is None:
    rep = self._rep                               # curver MappingClass objects directly
elif mod == -1:
    rep = self._rep.to_local_system()             # monodromy_rep.py:62 — integer homology
else:
    rep = self._rep.to_local_system_mod(mod)      # monodromy_rep.py:67 — mod-N homology
self._orbit = OrbitGraph(self._braid_gens, rep)   # orbit_graph.py:15
self._orbit.advance_until()                       # orbit_graph.py:139
return self._orbit
```

#### Visualization: `_draw(ax)` / `show_orbit()` / `save_orbit_png()`

`OrbitGraph.show()` is a stub (`orbit_graph.py:151`), so we implement our own.

```python
def _draw(self, ax):
    G = nx.convert_node_labels_to_integers(self._orbit.digraph)
    pos = nx.spring_layout(G, seed=42)            # reproducible layout
    # edges store generator as key= attribute (orbit_graph.py:76)
    edges = list(G.edges(keys=True))
    unique_keys = list(dict.fromkeys(k for _, _, k in edges))
    palette = plt.cm.tab10.colors
    color_map = {k: palette[i % 10] for i, k in enumerate(unique_keys)}
    edge_colors = [color_map[k] for _, _, k in edges]
    nx.draw_networkx(G, pos, ax=ax, arrows=True,
                     edge_color=edge_colors, node_size=300,
                     connectionstyle='arc3,rad=0.1')
    handles = [Line2D([0],[0], color=color_map[k], lw=2, label=f"σ_{k}")
               for k in unique_keys]
    ax.legend(handles=handles)

def show_orbit(self, figsize=(10, 8), layout='spring'):
    fig, ax = plt.subplots(figsize=figsize)
    self._draw(ax)
    plt.tight_layout()
    plt.show()

def save_orbit_png(self, filename, dpi=150):
    fig, ax = plt.subplots(figsize=(10, 8))
    self._draw(ax)
    plt.tight_layout()
    fig.savefig(filename, dpi=dpi)
    plt.close(fig)
    print(f"Saved to {filename}")
```

**Edge key detail**: `OrbitGraph.advance_node_gen` calls `add_edge(n, m, key=str(g))` (`orbit_graph.py:76`). After `convert_node_labels_to_integers`, each edge tuple is `(u, v, key_string)`.

#### `gap_generators()` / `gap_generators_file()`

Computes **stabilizer generators** of the orbit as **braid group words**, matching the approach in `LefschetzDisk.ipynb` (`build_stab_gens` / `gens_from_lcosets`). Requires `build_orbit(mod)` to have been called first — the `mod` is used only to make the orbit finite; the generators themselves are elements of B_m (no modular reduction applied to them).

**Algorithm** (port of the notebook's coset machinery to NetworkX + solomon):

```
1. BFS spanning tree from orbit.start in orbit.digraph
   → coset_rep[v]: FreeGrpElement word (path from start to v via braid gens)

2. For every directed edge (u → v, key=str(g)) in orbit.digraph:
   Schreier generator word = coset_rep[u] * g * coset_rep[v]^(-1)
   (This is the Schreier lemma: each edge gives one stabilizer generator candidate)

3. Discard identity words; deduplicate

4. Convert each FreeGrpElement Tietze list to a GAP braid word string:
   Tietze [1, 2, -1] → "s[1]*s[2]*s[1]^(-1)"

5. Emit GAP script defining the free group and subgroup
```

**Implementation sketch:**
```python
def gap_generators(self):
    assert self._orbit is not None, "Call build_orbit(mod) first"
    G = self._orbit.digraph
    start = self._orbit.start
    gen_map = {str(g): g for g in self._braid_gens}   # str → FreeGrpElement
    identity_word = self._braid_gens[0] ** 0

    # Step 1: BFS coset representatives
    coset_rep = {start: identity_word}
    queue = deque([start])
    while queue:
        u = queue.popleft()
        for v, edge_dict in G[u].items():
            if v not in coset_rep:
                key = next(iter(edge_dict))
                coset_rep[v] = coset_rep[u] * gen_map[key]
                queue.append(v)

    # Step 2: Schreier generator words
    stab_words = []
    for u, v, key in G.edges(keys=True):
        g = gen_map[key]
        word = coset_rep[u] * g * (coset_rep[v] ** -1)   # FreeGrpElement
        if word.Tietze():                                  # skip identity (empty Tietze)
            stab_words.append(word)

    # Step 3: deduplicate
    seen = []
    for w in stab_words:
        if not any(w == s for s in seen):
            seen.append(w)

    # Step 4: Tietze → GAP string
    def tietze_to_gap(tietze):
        if not tietze:
            return "Identity(F)"
        terms = []
        for i in tietze:
            if i > 0:
                terms.append(f"s[{i}]")
            else:
                terms.append(f"s[{-i}]^(-1)")
        return "*".join(terms)

    # Step 5: GAP script
    m = len(self._braid_gens)   # number of braid generators = m-1 for B_m
    lines = [tietze_to_gap(w.Tietze()) for w in seen]
    return (
        f"# Stabilizer generators in B_{m+1}"
        f"  (genus {self.genus}, {len(self._twists)} twists)\n"
        f"# Orbit computed mod {self._orbit.start.homology_mod_rep}\n"  # -1 if full curver orbit

        f"F := FreeGroup({m});\n"
        f"s := GeneratorsOfGroup(F);\n"
        f"stab_gens := [\n  " + ",\n  ".join(lines) + "\n];\n"
        f"H := Subgroup(F, stab_gens);\n"
    )

def gap_generators_file(self, filename):
    code = self.gap_generators()
    with open(filename, 'w') as f:
        f.write(code)
    print(f"GAP stabilizer generators written to {filename}")
```

**Correspondence to `LefschetzDisk.ipynb`:**
- BFS coset_rep = `build_lcosets(Gr, gens)` (cell 41)
- `coset_rep[u] * g * coset_rep[v]^(-1)` = Schreier generator from `gens_one_rcoset` (cell 41)
- Output is braid words, not matrices — the notebook's braid group elements are reproduced here as `FreeGrpElement` Tietze lists formatted for GAP
- The notebook's `braid_in_img` predicate (which further filters generators lying in the image of the monodromy map) is omitted — all stabilizer generators are exported; Ethan can apply additional filters in GAP

---

### `workspace.ipynb` — Cell Outline

| # | Content |
|---|---------|
| 1 | Imports: `from ethan import Fibration`, matplotlib inline setup |
| 2 | **Quick start**: 3-twist alternating genus-1 fibration (matches `alternating_genus_1(3)` from `examples.py`) |
| 3 | `build_orbit(mod=5)`, `show_orbit()` |
| 4 | `save_orbit_png('orbit.png')` |
| 5 | `gap_generators_file('gens.g')` — stabilizer generators as braid words via Schreier lemma (requires `build_orbit(mod=...)` first) |
| 6 | **Custom fibration**: user-chosen genus, punctures, and curver curves |
| 7 | **Named curves**: `fib.define_curve('alpha', S('a_0') * S('b_0'))` |
| 8 | **Loop over primes**: `for p in [2, 3, 5, 7]: fib.build_orbit(mod=p); ...` |

Each cell includes Markdown explaining what the cell does in plain English.

---

## Solomon Files Used

| File | How Used |
|------|----------|
| `src/solomon/monodromy_rep.py` | `MonodromyRep`, `to_local_system_mod`, `homology_matrices` |
| `src/solomon/orbit_graph.py` | `OrbitGraph`, `advance_until` |
| `src/solomon/surface_group.py` | `SurfaceGroup`, `FreeGrp` |
| `src/solomon/examples.py` | Pattern reference for domain group construction |

---

## Known Limitations (to document in notebook)

- `MonodromyRep.__call__` only supports modular inverse for 2×2 matrices (`monodromy_rep.py:86`). Mod-orbit computation for genus > 1 may fail if the braid action requires inverses mid-twist. Flag this if Ethan hits it.
- `FreeGrp(n-1).gens`: n Dehn twists → n−1 braid generators. This is the correct count (per `examples.py` pattern), but means you need at least 2 twists to have a non-trivial orbit.
- Custom mapping classes must live on `fib.surface` (the same curver surface). You cannot mix curves from differently-loaded surfaces.

---

## Verification Checklist

1. `pip install -e .` from solomon root, `cd <dir>`, launch JupyterLab
2. Run Example 1 → orbit size matches `alternating_genus_1(3)` from `examples.py`
3. `save_orbit_png('test.png')` → PNG file created on disk
4. `build_orbit(mod=5)` then `gap_generators_file('test.g')` → `.g` file with braid words; load in GAP and check `Index(F, H)`
5. Define custom curve via `S('a_0') * S('b_0')`, build orbit — no errors
6. Loop over `[2, 3, 5, 7]` primes — confirm different orbit sizes per prime
