"""Utility functions for sparse-grid basis and index operations."""

from __future__ import annotations


def cross(*args: list[list[int]]) -> list[list[int]]:
    """Compute pairwise cross-product (concatenation) of two index lists.

    Parameters
    ----------
    *args
        Two index collections to combine.

    Returns
    -------
    list[list[int]]
        Cross-product of index lists.
    """
    ans: list[list[int]] = []
    for arg in args[0]:
        for arg2 in args[1]:
            ans.append([*arg, *arg2])
    return ans


def eval_basis_1d(
    x: float, basis: list[int], interval: tuple[float, float] | None = None
) -> float:
    """Evaluate a one-dimensional hat basis function.

    Parameters
    ----------
    x
        Evaluation coordinate.
    basis
        Basis index as `[level, position]`.
    interval
        Optional physical interval. If omitted, `[0, 1]` is used.

    Returns
    -------
    float
        Value of 1-D hat basis at `x`.
    """
    if interval is None:
        return 1.0 - abs(x * 2 ** basis[0] - basis[1])
    pos = (x - interval[0]) / (interval[1] - interval[0])
    return 1.0 - abs(pos * 2 ** basis[0] - basis[1])
