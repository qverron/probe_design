"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from . import autocomplete, ifpd2
from . import db, extract_kmers
from . import query, query2

import logging
from rich.logging import RichHandler  # type: ignore

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler(markup=True, rich_tracebacks=True)],
)

__all__ = [
    "autocomplete",
    "ifpd2",
    "db",
    "extract_kmers",
    "query",
    "query2",
]
