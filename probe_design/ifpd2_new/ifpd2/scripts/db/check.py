"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import click  # type: ignore
from ...const import CONTEXT_SETTINGS
from ...database import DataBase
from ...walker2 import ChromosomeWalker
import logging
from tqdm import tqdm  # type: ignore


@click.command(
    name="check",
    context_settings=CONTEXT_SETTINGS,
    help="""Check integrity of INPUT database.""",
)
@click.argument("input_paths", metavar="INPUT", nargs=1, type=click.Path(exists=True))
def main(input_paths: str) -> None:
    DB = DataBase(input_paths)
    for chromosome in DB.chromosome_list:
        walker = ChromosomeWalker(DB, chromosome)
        previous_position = -1
        for record in tqdm(
            walker.buffer(),
            desc=f"Checking sorting '{chromosome.decode()}'",
            total=DB.chromosome_recordnos[chromosome],
        ):
            if record["start"] <= previous_position:
                raise AssertionError
            previous_position = record["start"]
    logging.info("That's all! :smiley:")
