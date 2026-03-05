# sparse_grid

**A Python Sparse Grid Package**

[![Tests](https://github.com/eggzec/sparse_grid/actions/workflows/code_test.yml/badge.svg)](https://github.com/eggzec/sparse_grid/actions/workflows/code_test.yml)
[![Documentation](https://github.com/eggzec/sparse_grid/actions/workflows/docs_build.yml/badge.svg)](https://github.com/eggzec/sparse_grid/actions/workflows/docs_build.yml)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

[![codecov](https://codecov.io/gh/eggzec/sparse-grid/branch/master/graph/badge.svg)](https://codecov.io/gh/eggzec/sparse-grid)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=eggzec_sparse_grid&metric=alert_status)](https://sonarcloud.io/project/overview?id=eggzec_sparse_grid)
[![License: BSD-3](https://img.shields.io/badge/License-BSD--3-blue.svg)](LICENSE)

[![PyPI Downloads](https://img.shields.io/pypi/dm/sparse-grid.svg?label=PyPI%20downloads)](https://pypi.org/project/sparse-grid/)
[![Python versions](https://img.shields.io/pypi/pyversions/sparse-grid.svg)](https://pypi.org/project/sparse-grid/)

`sparse_grid` is a pure-Python implementation of regular
[sparse grids](https://en.wikipedia.org/wiki/Sparse_grid)
over box domains. It provides hierarchical index generation, nodal-to-hierarchical
coefficient conversion, and fast function evaluation using the hat basis:

$$f_n(\mathbf{x}) = \sum_{\lvert\mathbf{l}\rvert_1 \le n+d-1} \sum_{\mathbf{p}} \alpha_{\mathbf{l},\mathbf{p}}\, \prod_{i=1}^{d} \phi_{l_i, p_i}(x_i)$$

## Quick example

```python
from sparse_grid import SparseGrid

sg = SparseGrid(dim=2, level=3)
sg.generate_points()

for index in sg.indices:
    pos = sg.g_p[tuple(index)].pos
    sg.g_p[tuple(index)].fv = (
        4.0 * pos[0] * (1.0 - pos[0]) * 4.0 * pos[1] * (1.0 - pos[1])
    )

sg.nodal_2_hier()
print(sg.eval_funct([0.25, 0.75]))
```

## Installation

```bash
pip install sparse_grid
```

Requires Python 3.10+. No external runtime dependencies. See the
[full installation guide](https://eggzec.github.io/sparse_grid/installation/) for
uv, poetry, and source builds.

## Documentation

- [Theory](https://eggzec.github.io/sparse_grid/theory/) — mathematical background, hierarchical basis, algorithms
- [Quickstart](https://eggzec.github.io/sparse_grid/quickstart/) — runnable examples
- [API Reference](https://eggzec.github.io/sparse_grid/api/) — class and function signatures
- [References](https://eggzec.github.io/sparse_grid/references/) — literature citations

## License

BSD-3-Clause — see [LICENSE](LICENSE).
