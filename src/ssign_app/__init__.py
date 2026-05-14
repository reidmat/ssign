"""ssign — Secretion-system Identification for Gram Negatives."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("ssign")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"
