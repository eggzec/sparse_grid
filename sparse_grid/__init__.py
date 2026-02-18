from importlib.metadata import PackageNotFoundError, version

from .grid import SparseGrid


try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = ["SparseGrid"]
