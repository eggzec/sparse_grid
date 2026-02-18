import unittest

from sparse_grid import SparseGrid


class TestFuncTest(unittest.TestCase):
    """simple test of sparse grid for sparse grid in 3d of level 3"""

    def test_sg_no_bound_3d(self) -> None:
        #
        #  SG is a sparse grid of dimension 3 and level 3.
        #  Create sg.indices which stores level and position for each point.
        #
        sg = SparseGrid(3, 3)
        #
        #  Determine sg.gP with the coordinates of the points
        #  associated with the sparse grid index set.
        #
        sg.generate_points()
        #
        #  Print the points in the grid.
        #
        print("")
        print("Coordinates of points in 3D sparse grid of level 3.")
        print("")
        for i in range(len(sg.indices)):
            sg.g_p[tuple(sg.indices[i])].print_point()
        #
        #  Did we compute the right number of grid points?
        #
        self.assertEqual(len(sg.indices), 31)
        #
        #  Evaluate 4x(1-x)*4y(1-y)*4z(1-z) at each grid point.
        #
        for i in range(len(sg.indices)):
            total = 1.0
            pos = sg.g_p[tuple(sg.indices[i])].pos
            for j in range(len(pos)):
                total *= 4.0 * pos[j] * (1.0 - pos[j])
            sg.g_p[tuple(sg.indices[i])].fv = total
        #
        #  Convert the sparse grid from nodal to hierarchical values.
        #
        sg.nodal_2_hier()
        #
        #  Does the evaluation of the sparse grid function in
        #  hierarchical values give the correct value gv?
        #
        for i in range(len(sg.indices)):
            self.assertEqual(
                sg.g_p[tuple(sg.indices[i])].fv,
                sg.eval_funct(sg.g_p[tuple(sg.indices[i])].pos),
            )

    def test_sg_no_bound_2d(self) -> None:
        #
        #  SG is a sparse grid of dimension 2 and level 3.
        #  Create sg.indices which stores level and position for each point.
        #
        sg = SparseGrid(2, 3)
        #
        #  Determine sg.gP with the coordinates of the points
        #  associated with the sparse grid index set.
        #
        sg.generate_points()
        #
        #  Print the points in the grid.
        #
        print("")
        print("Coordinates of points in 2D sparse grid of level 3.")
        print("")
        for i in range(len(sg.indices)):
            sg.g_p[tuple(sg.indices[i])].print_point()
        #
        #  Did we compute the right number of grid points?
        #
        self.assertEqual(len(sg.indices), 17)
        #
        #  Evaluate 4x(1-x)*4y(1-y) at each grid point.
        #
        for i in range(len(sg.indices)):
            total = 1.0
            pos = sg.g_p[tuple(sg.indices[i])].pos
            for j in range(len(pos)):
                total *= 4.0 * pos[j] * (1.0 - pos[j])
            sg.g_p[tuple(sg.indices[i])].fv = total
        #
        #  Convert to hierarchical values.
        #
        sg.nodal_2_hier()
        #
        #  Does the evaluation of sparse grid function in
        #  hierarchical values give the correct value gv?
        #
        for i in range(len(sg.indices)):
            self.assertEqual(
                sg.g_p[tuple(sg.indices[i])].fv,
                sg.eval_funct(sg.g_p[tuple(sg.indices[i])].pos),
            )
