# Theory

## 1) Full grids vs. sparse grids

A full tensor-product grid of level $n$ in $d$ dimensions has $O(2^{nd})$ points, making it impractical for $d > 3$. Sparse grids reduce this to

$$
\lvert\mathcal{H}_{n,d}\rvert = O\!\left(2^n \cdot n^{d-1}\right)
$$

points, while preserving accuracy for sufficiently smooth functions.

## 2) Hierarchical subspaces

The unit interval $[0,1]$ is decomposed into a hierarchy of subspaces indexed by level $l \ge 1$. At level $l$, the new grid points have positions

$$
x_{l,p} = \frac{p}{2^l}, \quad p \in \{1, 3, 5, \ldots, 2^l - 1\}.
$$

The number of new points at level $l$ is $2^{l-1}$.

In $d$ dimensions, the hierarchical subspace $W_{\mathbf{l}}$ for multi-index $\mathbf{l}=(l_1,\ldots,l_d)$ is the tensor product of one-dimensional subspaces:

$$
W_{\mathbf{l}} = W_{l_1} \otimes W_{l_2} \otimes \cdots \otimes W_{l_d}.
$$

## 3) Sparse grid admissibility

The regular sparse grid of level $n$ collects all subspaces satisfying:

$$
\lvert\mathbf{l}\rvert_1 = l_1 + l_2 + \cdots + l_d \le n + d - 1.
$$

This condition balances resolution across dimensions, discarding high-level subspaces that contribute little to smooth approximations.

## 4) Hat basis functions

The one-dimensional hat (piecewise-linear) basis function at level $l$, position $p$ is:

$$
\phi_{l,p}(x) = \max\!\left(0,\; 1 - \lvert 2^l x - p \rvert\right).
$$

It has support $\left[\frac{p-1}{2^l},\; \frac{p+1}{2^l}\right]$ and satisfies $\phi_{l,p}(x_{l,p}) = 1$.

The $d$-dimensional basis function is the product:

$$
\Phi_{\mathbf{l},\mathbf{p}}(\mathbf{x}) = \prod_{i=1}^{d} \phi_{l_i, p_i}(x_i).
$$

## 5) Hierarchical interpolation

Any function $f$ sampled at the sparse grid points can be represented as:

$$
f_n(\mathbf{x}) = \sum_{\lvert\mathbf{l}\rvert_1 \le n+d-1}\;\sum_{\mathbf{p}} \alpha_{\mathbf{l},\mathbf{p}}\;\Phi_{\mathbf{l},\mathbf{p}}(\mathbf{x}),
$$

where $\alpha_{\mathbf{l},\mathbf{p}}$ are the hierarchical surplus coefficients.

## 6) Nodal-to-hierarchical transform

Given nodal values $f(\mathbf{x}_{\mathbf{l},\mathbf{p}})$, the hierarchical surpluses are computed by the lifting algorithm. For each dimension $i$ and each level $l_i$ (processed from finest to coarsest), the one-dimensional transform subtracts the half-sum of neighboring coarser-level values:

$$
\alpha_{\mathbf{l},\mathbf{p}} = f_{\mathbf{l},\mathbf{p}} - \frac{1}{2}\!\left(f_{\mathbf{l}^-,\mathbf{p}^-} + f_{\mathbf{l}^-,\mathbf{p}^+}\right),
$$

where $(\mathbf{l}^-, \mathbf{p}^-)$ and $(\mathbf{l}^-, \mathbf{p}^+)$ are the left and right parent indices obtained by reducing the level and halving the position index. Boundary cases (where a parent falls outside the domain) use only the available neighbor.

## 7) Evaluation

To evaluate $f_n(\mathbf{x})$:

1. For each dimension $i$ and level $l_i = 1, \ldots, n$, determine the active basis index

$$
p_i = 2\left\lceil x_i \cdot 2^{l_i - 1}\right\rceil - 1
$$

and evaluate the corresponding hat function $\phi_{l_i, p_i}(x_i)$.

2. Iterate over all admissible hierarchical subspaces $\lvert\mathbf{l}\rvert_1 \le n + d - 1$.

3. For each subspace, multiply the per-dimension basis values and accumulate the weighted hierarchical coefficient:

$$
f_n(\mathbf{x}) \mathrel{+}= \alpha_{\mathbf{l},\mathbf{p}} \prod_{i=1}^{d} \phi_{l_i, p_i}(x_i).
$$

## 8) Domain scaling

For a general box domain $[a_i, b_i]$ instead of $[0,1]$, coordinates are mapped by

$$
\hat{x}_i = \frac{x_i - a_i}{b_i - a_i}
$$

before basis evaluation, and point positions are computed as

$$
x_i = a_i + (b_i - a_i)\,\frac{p_i}{2^{l_i}}.
$$
