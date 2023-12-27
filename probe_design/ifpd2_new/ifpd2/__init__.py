"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from . import asserts, dataclasses, fasta, io
from . import walker, walker2
from . import chromosome, database, region
from . import oligo, probe, probe_set
from . import loggingg

from importlib.metadata import version, PackageNotFoundError
from typing import List

__version__ = "unknown"

__all__ = [
    "asserts",
    "dataclasses",
    "fasta",
    "io",
    "walker",
    "walker2",
    "chromosome",
    "database",
    "region",
    "oligo",
    "probe",
    "probe_set",
    "loggingg",
]
__path__: List[str]

try:
    __version__ = version(__name__)
    __all__.append("__version__")
except PackageNotFoundError:
    pass

