# Getting Started

## Installation

```bash
pip install sparse_grid
```

## Create a sparse grid

```python
from sparse_grid import SparseGrid

sg = SparseGrid(dim=3, level=3)
sg.generate_points()
print(len(sg.indices))
```

## Assign nodal values

Populate nodal values (`fv`) at each sparse-grid point.

```python
for index in sg.indices:
    pos = sg.g_p[tuple(index)].pos
    value = 1.0
    for coord in pos:
        value *= 4.0 * coord * (1.0 - coord)
    sg.g_p[tuple(index)].fv = value
```

## Convert to hierarchical values

```python
sg.nodal_2_hier()
```

## Evaluate the sparse-grid interpolant

```python
x = [0.2, 0.4, 0.8]
y = sg.eval_funct(x)
print(y)
```
