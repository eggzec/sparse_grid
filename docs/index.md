# sparse_grid

`sparse_grid` is a compact Python implementation of regular sparse grids over
box domains.

It provides:

- sparse-grid index and point generation,
- nodal-to-hierarchical coefficient conversion,
- fast function evaluation using hierarchical basis functions.

## Documentation map

- **Getting Started**: installation and first working example.
- **User Guide**: core workflow and key data structures.
- **API Reference**: auto-generated module docs via `mkdocstrings`.

## Quick example

```python
from sparse_grid import SparseGrid

sg = SparseGrid(dim=2, level=3)
sg.generate_points()

for index in sg.indices:
    pos = sg.g_p[tuple(index)].pos
    sg.g_p[tuple(index)].fv = (
        4.0 * pos[0] * (1.0 - pos[0]) * 4.0 * pos[1] * (1.0 - pos[1])
    )

sg.nodal_2_hier()
value = sg.eval_funct([0.25, 0.75])
print(value)
```
