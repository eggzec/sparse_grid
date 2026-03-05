# API Reference

All public symbols are exported from the `sparse_grid` package.

## Exported symbols

- `SparseGrid`
- `GridPoint`
- `cross`
- `eval_basis_1d`

---

## `SparseGrid`

### `SparseGrid(dim=1, level=1)`

Create a regular sparse grid.

**Attributes**

- `dim` (`int`): number of dimensions.
- `level` (`int`): sparse-grid level.
- `indices` (`list[list[int]]`): generated multi-indices, each entry `[l_1, p_1, ..., l_d, p_d]`.
- `g_p` (`dict[tuple[int, ...], GridPoint]`): map from index tuple to `GridPoint`.
- `domain` (`tuple[tuple[float, float], ...]`): per-dimension bounds, default $[0, 1]^d$.

### `SparseGrid.generate_points()`

Generate all sparse-grid indices satisfying

$$
\lvert\mathbf{l}\rvert_1 \le n + d - 1
$$

and create the corresponding `GridPoint` objects stored in `g_p`.

### `SparseGrid.nodal_2_hier()`

Convert all nodal values (`fv`) to hierarchical surplus values (`hv`) using the dimension-wise lifting algorithm:

$$
\alpha_{\mathbf{l},\mathbf{p}} = f_{\mathbf{l},\mathbf{p}} - \frac{1}{2}\!\left(f_{\mathbf{l}^-,\mathbf{p}^-} + f_{\mathbf{l}^-,\mathbf{p}^+}\right).
$$

Must be called after all `fv` values are set.

### `SparseGrid.eval_funct(x) -> float`

Evaluate the sparse-grid interpolant at point `x`.

**Parameters**

- `x` (`list[float]`): evaluation point in physical coordinates, length must equal `dim`.

**Returns**

The interpolated value

$$
f_n(\mathbf{x}) = \sum_{\lvert\mathbf{l}\rvert_1 \le n+d-1}\;\sum_{\mathbf{p}} \alpha_{\mathbf{l},\mathbf{p}}\;\prod_{i=1}^{d}\phi_{l_i,p_i}(x_i).
$$

Must be called after `nodal_2_hier()`.

### `SparseGrid.print_grid()`

Print the current hierarchical subspace index to stdout.

### `SparseGrid.loop_hier_spaces()`

Iterate over all hierarchical subspaces satisfying the admissibility condition and invoke the current `action` callback for each.

### `SparseGrid.generate_points_rec(dim, level, cur_level=None) -> list[list[int]]`

Recursively generate multi-indices for all admissible hierarchical subspaces.

**Parameters**

- `dim` (`int`): remaining dimensions.
- `level` (`int`): remaining level budget.
- `cur_level` (`int | None`): current recursion level (default: `1`).

**Returns**

List of interleaved multi-indices `[l_1, p_1, ..., l_d, p_d]`.

### `SparseGrid.nodal_2_hier_1d(node, i, j, dim)`

Apply the one-dimensional nodal-to-hierarchical transform along dimension `dim` at level `i`, position `j`.

**Parameters**

- `node` (`list[int]`): fixed index components from other dimensions.
- `i` (`int`): level index.
- `j` (`int`): position index.
- `dim` (`int`): dimension being transformed.

---

## `GridPoint`

### `GridPoint(index=None, domain=None)`

Representation of a single sparse-grid point.

**Attributes**

- `pos` (`list[float]`): physical coordinates.
- `fv` (`float`): nodal function value (default `0.0`).
- `hv` (`float`): hierarchical surplus value (default `0.0`).

### `GridPoint.point_position(index, domain=None) -> list[float]`

Convert a multi-index to physical coordinates.

**Parameters**

- `index` (`list[int]`): interleaved multi-index `[l_1, p_1, ..., l_d, p_d]`.
- `domain` (`tuple[tuple[float, float], ...] | None`): per-dimension bounds.

Without a domain the mapping is:

$$
x_i = \frac{p_i}{2^{l_i}}.
$$

With domain bounds $[a_i, b_i]$:

$$
x_i = a_i + (b_i - a_i)\,\frac{p_i}{2^{l_i}}.
$$

### `GridPoint.print_point()`

Print point coordinates as a tab-separated line to stdout.

---

## Utility functions

### `cross(*args) -> list[list[int]]`

Compute the pairwise cross-product (concatenation) of two index lists. Used internally by `generate_points_rec` to combine per-dimension index sets.

**Parameters**

- `*args`: two index collections `list[list[int]]`.

**Returns**

Cross-product of the two lists.

### `eval_basis_1d(x, basis, interval=None) -> float`

Evaluate a one-dimensional hat basis function:

$$
\phi_{l,p}(x) = \max\!\left(0,\; 1 - \lvert 2^l x - p \rvert\right).
$$

**Parameters**

- `x` (`float`): evaluation coordinate.
- `basis` (`list[int]`): basis index `[level, position]`.
- `interval` (`tuple[float, float] | None`): physical interval (default: $[0, 1]$).

**Returns**

Value of the hat basis at `x`.

## References

See [References](references.md) for full citations.
