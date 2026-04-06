# Solomon's Desk

A JupyterLab workspace for computing orbit graphs and stabilizer generators of Lefschetz fibrations. Built on top of the **solomon** library.

Open `workspace.ipynb` to get started.

---

## Setup

From the **solomon root directory**, launch JupyterLab using `uv`:

```
./run_jupyter
```

Then navigate to `solomons_desk/workspace.ipynb` in the file browser and open it.

All dependencies (curver, networkx, matplotlib, …) are managed automatically by `uv` — no manual install step needed.

---

## Quick Start

```python
from ethan import Fibration

# 1. Create a fibration: genus-1 fiber surface with 1 puncture
fib = Fibration(genus=1, punctures=1)

# 2. Add Dehn twist generators by curve name
fib.add_twist('a_0')
fib.add_twist('b_0')
fib.add_twist('a_0')

# 3. Compute the orbit graph mod 5 (makes it finite)
fib.build_orbit(mod=5)

# 4. View the graph
fib.show_orbit()

# 5. Export stabilizer generators to GAP
fib.gap_generators_file('my_stabilizer.g')
```

---

## `Fibration` — Method Reference

### Creating a fibration

```python
fib = Fibration(genus, punctures)
```

- `genus` — genus of the fiber surface (e.g. `1` for a torus)
- `punctures` — number of punctures on the fiber (use `1` for one boundary component)

The underlying curver surface is available as `fib.surface`.

---

### Adding Dehn twists

```python
fib.add_twist('a_0')          # standard curver curve name
fib.add_twist('b_0')
fib.add_twist(fib.surface('a_0') * fib.surface('b_0'))  # composition directly
```

You need **at least 2 twists** before computing an orbit.

Available standard curve names on a genus-*g* surface:

| Name | Description |
|------|-------------|
| `a_0`, `a_1`, …, `a_{g-1}` | Handle *a*-curves |
| `b_0`, `b_1`, …, `b_{g-1}` | Handle *b*-curves |
| `c_0`, `c_1`, …, `c_{g-2}` | Chain curves between handles |

---

### Naming custom mapping classes

```python
fib.define_curve('name', curve)
```

- `curve` can be a string (`'a_0'`) or a curver MappingClass object
- Once defined, use the name with `add_twist`

```python
S = fib.surface
fib.define_curve('conj_a', S('a_0') * S('b_0') * S('a_0') ** -1)
fib.add_twist('conj_a')
```

---

### Building the orbit graph

```python
fib.build_orbit(mod=5)    # finite orbit mod 5  ← most common
fib.build_orbit(mod=-1)   # integer homology matrices (may not terminate)
fib.build_orbit()         # full curver representation (may not terminate)
```

| `mod` value | What it does |
|-------------|--------------|
| positive int (e.g. `5`) | Reduces representation mod N → finite orbit |
| `-1` | Converts to integer homology matrices |
| `None` (default) | Uses curver mapping class objects directly |

After building, check the size:

```python
fib.orbit_size()   # → number of nodes
```

---

### Viewing the orbit graph

```python
fib.show_orbit()                         # display inline in Jupyter
fib.show_orbit(figsize=(14, 10))         # larger figure
fib.save_orbit_png('orbit.png')          # save to file
fib.save_orbit_png('orbit.png', dpi=200) # higher resolution
```

Edges are **colour-coded by generator** (σ₁, σ₂, …). A legend is shown.

---

### Exporting to GAP

```python
fib.gap_generators()               # returns a string
fib.gap_generators_file('stab.g')  # writes to a file
```

**Important:** you must call `build_orbit(mod=N)` *before* this. The mod-N reduction is what makes the orbit finite and gives a meaningful stabilizer.

The output is a GAP script defining the stabilizer subgroup of the braid group:

```gap
F := FreeGroup(2);
s := GeneratorsOfGroup(F);
stab_gens := [
  s[1]*s[2]*s[1]^(-1),
  s[2]*s[1]*s[2]^(-1)
];
H := Subgroup(F, stab_gens);
```

Here `s[i]` is the standard Artin generator σᵢ. The generators are computed using the **Schreier lemma** applied to the orbit graph: for every edge `u →σᵢ→ v`, the element `t_u · σᵢ · t_v⁻¹` (where `t_w` is a spanning-tree path from the base vertex to `w`) lies in the stabilizer.

---

## Looping over primes

```python
for p in [2, 3, 5, 7, 11]:
    fib.build_orbit(mod=p)
    print(f"mod {p}: {fib.orbit_size()} nodes")
    fib.save_orbit_png(f'orbit_mod{p}.png')
    fib.gap_generators_file(f'stab_mod{p}.g')
```

---

## Known limitations

- **Genus > 1 with mod orbit:** the modular matrix inverse is only implemented for 2×2 matrices (genus-1 fiber). If you hit an error with higher genus, use `build_orbit(mod=-1)` (integer matrices) or `build_orbit()` (full curver).
- **Infinite orbits:** `build_orbit()` and `build_orbit(mod=-1)` may not terminate for generic fibrations. Always try `mod=N` first.
- **Minimum 2 twists:** the braid group on *n* marked points needs *n* ≥ 2 twists (giving σ₁, …, σₙ₋₁).
- **Custom curves must share the same surface:** all mapping classes passed to `add_twist` must come from `fib.surface`. You cannot mix curves from different surface objects.
