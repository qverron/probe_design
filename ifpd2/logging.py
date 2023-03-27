"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import logging
from rich.logging import RichHandler  # type: ignore
from rich.console import Console  # type: ignore
import os
from typing import Optional


def add_log_file_handler(path: str, logger_name: Optional[str] = None) -> None:
    if os.path.isdir(path):
        raise AssertionError
    log_dir = os.path.dirname(path)
    if not (os.path.isdir(log_dir) or log_dir == ""):
        raise AssertionError
    fh = RichHandler(console=Console(file=open(path, mode="w+")), markup=True)
    fh.setLevel(logging.INFO)
    logging.getLogger(logger_name).addHandler(fh)
    logging.getLogger(logger_name).info(f"[green]Log to[/]: '{path}'")
