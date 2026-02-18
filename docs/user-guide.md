# User Guide

## Core workflow

1. Build a `SparseGrid(dim, level)`.
2. Generate points with `generate_points()`.
3. Set nodal values in each point's `fv`.
4. Convert nodal values using `nodal_2_hier()`.
5. Evaluate with `eval_funct(x)`.

## Main objects

- `SparseGrid`
  - Manages sparse-grid structure and transforms.
- `GridPoint`
  - Stores point coordinates (`pos`), nodal value (`fv`), and hierarchical value (`hv`).

## Data layout

- `indices`: list of interleaved indices `[l_1, p_1, ..., l_d, p_d]`.
- `g_p`: dictionary keyed by index tuple, value is `GridPoint`.
- `domain`: tuple of intervals per dimension, default `[0, 1]^d`.

## Typical usage notes

- Always call `generate_points()` before filling values.
- Always call `nodal_2_hier()` before `eval_funct()`.
- Use consistent dimension between `x` and the grid `dim`.
