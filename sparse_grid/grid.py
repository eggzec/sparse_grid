"""Sparse-grid container and algorithms."""

from __future__ import annotations

import copy
import math

from .point import GridPoint
from .utils import cross, eval_basis_1d


EVEN_DIVISOR = 2


class SparseGrid:
    """Regular sparse grid over a box domain.

    A grid is defined by:
    - `dim`: number of dimensions,
    - `level`: sparse-grid level,
    - `indices`: generated hierarchical indices,
    - `g_p`: map from index tuple to `GridPoint`.
    """

    def __init__(self, dim: int = 1, level: int = 1) -> None:
        self.dim = dim
        self.level = level
        self.g_p: dict[tuple[int, ...], GridPoint] = {}
        self.indices: list[list[int]] = []
        self.domain = ((0.0, 1.0),) * dim
        self.action = self.eval_action

    def print_grid(self) -> None:
        """Print the currently active hierarchical subspace index."""
        print(self.h_space)

    def eval_action(self) -> None:
        """Accumulate contribution of the current hierarchical subspace."""
        basis = copy.deepcopy(self.eval_per_dim[0][self.h_space[0] - 1][0])
        value = self.eval_per_dim[0][self.h_space[0] - 1][1]
        for i in range(1, self.dim):
            value *= self.eval_per_dim[i][self.h_space[i] - 1][1]
            basis += self.eval_per_dim[i][self.h_space[i] - 1][0]
        self.value += self.g_p[tuple(basis)].hv * value

    def eval_funct(self, x: list[float]) -> float:
        """Evaluate a sparse grid function.

        Hierarchical values have to be set.

        Parameters
        ----------
        x
            Evaluation point in physical coordinates.

        Returns
        -------
        float
            Evaluated sparse-grid function value at `x`.
        """
        self.value = 0.0
        self.eval_per_dim: list[list[list[list[int]] | float]] = []
        for i in range(self.dim):
            self.eval_per_dim.append([])
            for j in range(1, self.level + 1):
                pos = (x[i] - self.domain[i][0]) / (
                    self.domain[i][1] - self.domain[i][0]
                )
                basis = int(math.ceil(pos * 2 ** (j - 1)) * 2 - 1)
                if basis == -1:
                    basis = 1
                    self.eval_per_dim[i].append([[j, basis]])
                else:
                    self.eval_per_dim[i].append([[j, basis]])
                self.eval_per_dim[i][j - 1].append(
                    eval_basis_1d(
                        x[i], self.eval_per_dim[i][j - 1][0], self.domain[i]
                    )
                )
        self.action = self.eval_action
        self.loop_hier_spaces()
        return self.value

    def loop_hier_spaces(self) -> None:
        """Iterate over all hierarchical subspaces of the sparse grid."""
        for i in range(1, self.level + 1):
            self.h_space = [i]
            self.loop_hier_spaces_rec(self.dim - 1, self.level - (i - 1))

    def loop_hier_spaces_rec(self, dim: int, level: int) -> None:
        """Recursively traverse hierarchical subspaces.

        Parameters
        ----------
        dim
            Remaining dimensions to recurse through.
        level
            Remaining level budget for sparse-grid admissibility.
        """
        if dim > 1:
            for i in range(1, level + 1):
                self.h_space.append(i)
                self.loop_hier_spaces_rec(dim - 1, level - (i - 1))
                self.h_space.pop()
        else:
            for i in range(1, level + 1):
                self.h_space.append(i)
                self.action()
                self.h_space.pop()

    def generate_points(self) -> None:
        """Generate sparse-grid indices and corresponding `GridPoint` values."""
        self.indices = self.generate_points_rec(self.dim, self.level)
        for i in range(len(self.indices)):
            self.g_p[tuple(self.indices[i])] = GridPoint(
                self.indices[i], self.domain
            )

    def generate_points_rec(
        self, dim: int, level: int, cur_level: int | None = None
    ) -> list[list[int]]:
        """Run over all hierarchical subspaces and add all their indices.

        Returns
        -------
        list[list[int]]
            All generated sparse-grid multi-indices.
        """
        basis_cur = []
        if cur_level is None:
            cur_level = 1
        for i in range(1, 2 ** (cur_level) + 1, 2):
            basis_cur.append([cur_level, i])
        if dim == 1 and cur_level == level:
            return basis_cur
        if dim == 1:
            basis_cur += self.generate_points_rec(dim, level, cur_level + 1)
            return basis_cur
        if cur_level == level:
            return cross(
                basis_cur,
                self.generate_points_rec(dim - 1, level - cur_level + 1),
            )
        return cross(
            basis_cur, self.generate_points_rec(dim - 1, level - cur_level + 1)
        ) + self.generate_points_rec(dim, level, cur_level + 1)

    def nodal_2_hier_1d(
        self, node: list[int], i: int, j: int, dim: int
    ) -> None:
        """Apply one-dimensional nodal-to-hierarchical transform.

        Parameters
        ----------
        node
            Fixed index components from all dimensions except `dim`.
        i
            Level index in transformed dimension.
        j
            Position index in transformed dimension.
        dim
            Dimension currently being transformed.
        """
        left = [i - 1, j // 2]
        right = [i - 1, j // 2 + 1]
        while left[1] % EVEN_DIVISOR == 0 and left[0] > 0:
            left = [left[0] - 1, left[1] // EVEN_DIVISOR]
        while right[1] % EVEN_DIVISOR == 0 and right[0] > 0:
            right = [right[0] - 1, right[1] // EVEN_DIVISOR]
        if len(node) > EVEN_DIVISOR:
            pre_cur_dim = node[0 : 2 * dim]
            post_cur_dim = node[2 * dim : len(node) + 1]
            index = [*pre_cur_dim, i, j, *post_cur_dim]
            left = [*pre_cur_dim, *left, *post_cur_dim]
            right = [*pre_cur_dim, *right, *post_cur_dim]
        elif dim == 0:
            index = [i, j, *node]
            left += node
            right += node
        else:
            index = [*node, i, j]
            left = [*node, *left]
            right = [*node, *right]

        if left[2 * dim] == 0:
            if right[2 * dim] != 0:
                self.g_p[tuple(index)].hv -= 0.5 * self.g_p[tuple(right)].hv
        elif right[2 * dim] == 0:
            self.g_p[tuple(index)].hv -= 0.5 * self.g_p[tuple(left)].hv
        else:
            self.g_p[tuple(index)].hv -= 0.5 * (
                self.g_p[tuple(left)].hv + self.g_p[tuple(right)].hv
            )

    def nodal_2_hier(self) -> None:
        """Convert all nodal values in the grid to hierarchical values."""
        for i in range(len(self.indices)):
            self.g_p[tuple(self.indices[i])].hv = self.g_p[
                tuple(self.indices[i])
            ].fv
        for d in range(self.dim):
            for i in range(self.level, 0, -1):
                indices = self.generate_points_rec(
                    self.dim - 1, self.level - i + 1
                )
                for j in range(1, 2**i + 1, 2):
                    for k in range(len(indices)):
                        self.nodal_2_hier_1d(indices[k], i, j, d)
