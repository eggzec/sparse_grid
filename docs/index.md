# sparse_grid

A pure-Python implementation of regular sparse grids over box domains.

## Features

- Sparse-grid index and point generation via recursive Smolyak construction.
- Nodal-to-hierarchical coefficient conversion using the lifting algorithm.
- Fast function evaluation on hierarchical hat-basis representations.
- Arbitrary box domains $[a_1, b_1] \times \cdots \times [a_d, b_d]$.

## Mathematical model

A regular sparse grid of level $n$ in $d$ dimensions is the union of hierarchical subspaces satisfying the admissibility condition

$$
\mathcal{H}_{n,d} = \bigcup_{\lvert\mathbf{l}\rvert_1 \le n + d - 1} W_{\mathbf{l}},
$$

where $\mathbf{l} = (l_1, \ldots, l_d)$ is a level multi-index and $\lvert\mathbf{l}\rvert_1 = l_1 + \cdots + l_d$.

Each subspace $W_{\mathbf{l}}$ contains points indexed by odd positions $p_i \in \{1, 3, \ldots, 2^{l_i} - 1\}$. The physical coordinates on $[0, 1]^d$ are

$$
x_i = \frac{p_i}{2^{l_i}}, \quad i = 1, \ldots, d.
$$

The interpolant is expressed in the hierarchical hat basis

$$
f_n(\mathbf{x}) = \sum_{\lvert\mathbf{l}\rvert_1 \le n+d-1} \sum_{\mathbf{p}} \alpha_{\mathbf{l},\mathbf{p}}\, \prod_{i=1}^{d} \phi_{l_i, p_i}(x_i),
$$

where $\alpha_{\mathbf{l},\mathbf{p}}$ are hierarchical surplus coefficients and $\phi_{l,p}$ is the one-dimensional hat basis function.

## Public API

- `SparseGrid` — grid container, point generation, transforms, and evaluation
- `GridPoint` — point storage with coordinates and function values
- `cross` — index cross-product utility
- `eval_basis_1d` — one-dimensional hat basis evaluation

## Documentation

- [Theory](theory.md) — mathematical background and algorithms
- [Installation](installation.md) — install from PyPI or Git
- [Quickstart](quickstart.md) — runnable examples
- [API Reference](api.md) — class and function signatures
- [References](references.md) — source literature
