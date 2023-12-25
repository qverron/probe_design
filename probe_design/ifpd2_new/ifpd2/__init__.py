"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from ifpd2 import asserts, dataclasses, fasta, io
from ifpd2 import walker, walker2
from ifpd2 import chromosome, database, region
from ifpd2 import oligo, probe, probe_set

from importlib.metadata import version, PackageNotFoundError
from typing import List

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    pass

__all__ = [
    "__version__",
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
]
__path__: List[str]
