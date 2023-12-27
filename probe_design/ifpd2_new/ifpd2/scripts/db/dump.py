"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import click  # type: ignore
from ... import const
from ...database import DataBase
from ...walker2 import ChromosomeWalker
from tqdm import tqdm  # type: ignore
from typing import List, Optional, Tuple


@click.command(
    name="dump",
    context_settings=const.CONTEXT_SETTINGS,
    help="""Dump whole INPUT database to the terminal, in tsv format.
    To dump a single feature, specify it with the --chrom option.
    To dump a region of a feature, also use the --start and --end options.""",
)
@click.argument("input_paths", metavar="INPUT", nargs=1, type=click.Path(exists=True))
@click.option("--chrom", type=click.STRING, help="Database chromosome to dump from.")
@click.option("--start", type=click.INT, help="Start location of the region to dump.")
@click.option("--end", type=click.INT, help="End location of the region to dump.")
def main(
    input_paths: str,
    chrom: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
) -> None:
    DB = DataBase(input_paths)
    chrom, start, end = check_region(DB, chrom, start, end)

    print("\t".join(const.database_columns))
    for chromosome in get_chromosome_list(DB, chrom):
        walker = ChromosomeWalker(DB, chromosome)
        for record in tqdm(
            walker.buffer(start, end),
            desc=f"dumping '{chromosome.decode()}'",
            total=DB.chromosome_recordnos[chromosome],
        ):
            if start > record["start"]:
                continue
            print(record.to_csv("\t"))


def check_region(
    DB: DataBase,
    chrom: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
) -> Tuple[Optional[str], int, int]:
    if chrom is None:
        if start is not None or end is not None:
            raise AssertionError(
                "cannot use --region-start or --region-end without --chrom"
            )
    elif start is not None:
        chrom_size = DB.chromosome_sizes_nt[chrom.encode()]
        if start >= chrom_size:
            raise AssertionError(f"{start} larger than chromosome size: {chrom_size}")
        if end is not None and start >= end:
            raise AssertionError("end location smaller than start")

    return (
        chrom,
        0 if start is None else start,
        -1 if end is None else end,
    )


def get_chromosome_list(DB: DataBase, chrom: Optional[str]) -> List[bytes]:
    chromosome_list = DB.chromosome_list
    if chrom is not None:
        if chrom.encode() not in chromosome_list:
            raise AssertionError(f"'{chrom}' not found")
        chromosome_list = [chrom.encode()]
    return chromosome_list
