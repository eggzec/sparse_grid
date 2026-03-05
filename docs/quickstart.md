# Quickstart

## Example 1: 2-D sparse grid interpolation

Interpolate the product function

$$
f(x, y) = 4x(1-x)\cdot 4y(1-y)
$$

on a level-3 sparse grid over $[0,1]^2$ with $17$ points.

```py
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
print(f"Number of grid points: {len(sg.indices)}")
print(f"f(0.25, 0.75) ≈ {value}")
```

Expected output:

```text
Number of grid points: 17
f(0.25, 0.75) ≈ 0.5625
```

## Example 2: 3-D sparse grid interpolation

Interpolate the product function

$$
f(x, y, z) = \prod_{i=1}^{3} 4\,x_i(1 - x_i)
$$

on a level-3 sparse grid over $[0,1]^3$ with $31$ points.

```py
from sparse_grid import SparseGrid

sg = SparseGrid(dim=3, level=3)
sg.generate_points()

for index in sg.indices:
    pos = sg.g_p[tuple(index)].pos
    total = 1.0
    for coord in pos:
        total *= 4.0 * coord * (1.0 - coord)
    sg.g_p[tuple(index)].fv = total

sg.nodal_2_hier()

# Evaluate at a test point and verify against the exact value
x = [0.2, 0.4, 0.8]
approx = sg.eval_funct(x)
exact = 1.0
for xi in x:
    exact *= 4.0 * xi * (1.0 - xi)
print(f"Number of grid points: {len(sg.indices)}")
print(f"Approximate: {approx}")
print(f"Exact:       {exact}")
```

Expected output:

```text
Number of grid points: 31
Approximate: 0.39321600000000013
Exact:       0.39321600000000013
```

## Example 3: evaluating at every grid point

Verify that the hierarchical interpolant exactly recovers all nodal values:

```py
from sparse_grid import SparseGrid

sg = SparseGrid(dim=2, level=3)
sg.generate_points()

for index in sg.indices:
    pos = sg.g_p[tuple(index)].pos
    sg.g_p[tuple(index)].fv = (
        4.0 * pos[0] * (1.0 - pos[0]) * 4.0 * pos[1] * (1.0 - pos[1])
    )

sg.nodal_2_hier()

max_err = 0.0
for index in sg.indices:
    pos = sg.g_p[tuple(index)].pos
    fv = sg.g_p[tuple(index)].fv
    approx = sg.eval_funct(pos)
    max_err = max(max_err, abs(fv - approx))

print(f"Max error at grid points: {max_err}")
```

Expected output:

```text
Max error at grid points: 0.0
```
