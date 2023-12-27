"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import click  # type: ignore
from ...const import CONTEXT_SETTINGS
from ...database import DataBase
import logging


@click.command(
    name="info", context_settings=CONTEXT_SETTINGS, help="Show INPUT database details."
)
@click.argument("input_paths", metavar="INPUT", nargs=1, type=click.Path(exists=True))
def main(input_paths: str) -> None:
    DB = DataBase(input_paths)
    DB.log_details()
    logging.info("")
    logging.info("That's all! :smiley:")
