"""Point model used by sparse-grid containers."""

from __future__ import annotations


class GridPoint:
    """Representation of a sparse-grid point.

    Stores the Cartesian position and both nodal (`fv`) and hierarchical
    (`hv`) function values.
    """

    def __init__(
        self,
        index: list[int] | None = None,
        domain: tuple[tuple[float, float], ...] | None = None,
    ) -> None:
        self.hv: float = 0.0
        self.fv: float = 0.0
        if index is None:
            self.pos: list[float] = []
        else:
            self.pos = self.point_position(index, domain)

    @staticmethod
    def point_position(
        index: list[int], domain: tuple[tuple[float, float], ...] | None = None
    ) -> list[float]:
        """Convert a level/position index into physical coordinates.

        Parameters
        ----------
        index
            Interleaved multi-index `[l_1, p_1, ..., l_d, p_d]`.
        domain
            Optional domain bounds per dimension. If omitted, `[0, 1]^d`
            is assumed.

        Returns
        -------
        list[float]
            Point coordinates in physical space.
        """
        coord = []
        if domain is None:
            for i in range(len(index) // 2):
                coord.append(index[2 * i + 1] / 2.0 ** index[2 * i])
        else:
            for i in range(len(index) // 2):
                coord.append(
                    (domain[i][1] - domain[i][0])
                    * index[2 * i + 1]
                    / 2.0 ** index[2 * i]
                    + domain[i][0]
                )
        return coord

    def print_point(self) -> None:
        """Print point coordinates as a tab-separated line."""
        if not self.pos:
            return
        out = ""
        for i in range(len(self.pos)):
            out += str(self.pos[i]) + "\t"
        print(out)
