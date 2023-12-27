"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import click  # type: ignore
from ..const import CONTEXT_SETTINGS
from ..dataclasses import Folder
from ..dataclasses import PositiveInteger
from ..dataclasses import QueryWindow
from ..loggingg import add_log_file_handler
from ..database import DataBase
from ..region import GenomicRegionBuilder
from ..walker2 import ChromosomeWalker
import logging
from os.path import isdir, isfile
from os import mkdir
from typing import Optional, Tuple


@click.command(
    name="query2",
    context_settings=CONTEXT_SETTINGS,
    help="""Design FISH probes on a database chromosome.
    When a --region is not specified, the whole chromosome is used.
    Use -M for the number of probes, and -N for the number of oligos per probe.
    Use -D to specify the minimum distance between consecutive oligos in a probe.

    Use --single to access single-probe design. Same as "-X 1 -w 1", this option
    overrides -X and -w, and ignores -W, if specified.""",
)
@click.argument("db_folder", type=click.Path(exists=True))
@click.argument("chromosome", type=click.STRING)
@click.argument("output_path", metavar="OUTPUT", type=click.Path(exists=False))
@click.option(
    "--region",
    type=click.INT,
    nargs=2,
    default=(0, -1),
    help="Space-separate region start and end location.",
)
@click.option("-M", "--probes", type=click.INT, help="Design M probes.")
@click.option(
    "-N",
    "--oligos",
    type=click.INT,
    default=48,
    help="N oligos per probe. Default: 48",
)
@click.option(
    "-D",
    "--dist",
    type=click.INT,
    default=2,
    help="Consecutive oligo distance. Default: 2",
)
@click.option("--single", is_flag=True, help="Easy access to single probe design.")
@click.option("-W", "--window-size", type=click.INT)
@click.option(
    "-w",
    "--window-shift",
    type=click.FLOAT,
    default=0.1,
    help="Shift as window fraction. Default: 0.1",
)
@click.option(
    "-R",
    "--focus-size",
    type=click.FLOAT,
    default=8e3,
    help="""\b
    Focus size in nt or window fraction.
    Default: 8000""",
)
@click.option(
    "-r",
    "--focus-step",
    type=click.FLOAT,
    default=1e3,
    help="""\b
    Focus step in nt or focus fraction.
    Default: 1000""",
)
@click.option(
    "-F",
    "--off-targets",
    type=click.INT,
    nargs=2,
    default=(0, 99),
    help="""Extremities of acceptable off-target range.
    Default: [0, 99]""",
)
@click.option(
    "-G",
    "--free-energy",
    type=click.FLOAT,
    nargs=2,
    default=(0.0, 0.5),
    help="""Extremities of acceptable 2nd structure dG range,
    in absolute kcal/mol or hybridization dG fraction.
    Default: [0.0, 0.5]""",
)
@click.option(
    "-o",
    "--oligo-score-step",
    type=click.FLOAT,
    default=0.1,
    help="Oligo score relaxation step. Default: 0.1",
)
@click.option(
    "-t",
    "--melting-half-width",
    type=click.FLOAT,
    default=10.0,
    help="Melting range half-width. Default: 10.0",
)
@click.option(
    "-P",
    "--probe-size",
    type=click.INT,
    default=1e4,
    help="Max probe size in nt. Default: 10000",
)
@click.option(
    "-H",
    "--gap-size",
    type=click.FLOAT,
    default=0.1,
    help="""\b
    Max gap size as probe size fraction.
    Default: 0.1""",
)
@click.option(
    "-I",
    "--oligo-intersection",
    type=click.FLOAT,
    default=0.5,
    help="""Probe intersection threshold as shared oligo fraction.
    Default: .5""",
)
@click.option(
    "-k",
    "--oligo-length",
    type=click.INT,
)
@click.option(
    "--threads",
    type=click.INT,
    default=1,
    help="Threads for parallelization. Default: 1",
)
def main(
    db_folder: str,
    chromosome: str,
    output_path: str,
    region: Tuple[int, int],
    probes: Optional[int],
    oligos: int,
    dist: int,
    single: bool,
    window_size: Optional[int],
    window_shift: Optional[float],
    focus_size: float,
    focus_step: float,
    off_targets: Tuple[int, int],
    free_energy: Tuple[float, float],
    oligo_score_step: float,
    melting_half_width: float,
    probe_size: int,
    gap_size: float,
    oligo_intersection: float,
    oligo_length: Optional[int],
    threads: int,
) -> None:
    settings = QuerySettings(db_folder, output_path)

    probes, window_size, window_shift = check_input(
        single, probes, window_size, window_shift
    )

    if not isdir(settings.output_path):
        mkdir(settings.output_path)
    add_log_file_handler(f"{settings.output_path}/ifpd2-main.log")

    DB = DataBase(settings.db_folder_path)

    RB = GenomicRegionBuilder(DB.get_chromosome(chromosome.encode()))
    if probes is not None:
        region_set_list = RB.build_by_number(PositiveInteger(probes).n)
    else:
        window_size, window_shift = QueryWindow(window_size, window_shift).astuple()
        if window_size is None or window_shift is None:
            raise AssertionError
        region_set_list = RB.build_by_size(window_size, window_shift)

    walker = ChromosomeWalker(DB, chromosome.encode())
    walker.walk_multiple_regions(region_set_list)

    # print((DB, RB, walker, region_set_list))

    logging.info("Done. :thumbs_up: :smiley:")
    logging.shutdown()


def check_input(
    single: bool,
    probes: Optional[int],
    window_size: Optional[int],
    window_shift: Optional[float],
) -> Tuple[Optional[int], Optional[int], Optional[float]]:
    if single:
        probes = 1
        window_shift = None
        window_size = None

    if probes is None and window_size is None:
        raise AssertionError
    if probes is not None and window_size is not None:
        logging.warning("cannot combine -X and -W. Using -X.")
        window_size = None
        window_shift = None
    if window_size is not None and window_shift is None:
        window_shift = 1

    return (probes, window_size, window_shift)


class QuerySettings:
    _db_folder: Folder
    _output_path: Folder

    def __init__(self, db_folder: str, output_path: str):
        super(QuerySettings, self).__init__()
        self._db_folder = Folder(db_folder, True)
        self._output_path = Folder(output_path, False)

    @property
    def db_folder_path(self) -> str:
        return self._db_folder.path

    @property
    def output_path(self) -> str:
        return self._output_path.path

    @property
    def reusable(self) -> bool:
        return False


def assert_reusable(output_path: str):
    assert_msg = " ".join(
        [
            "output path should NOT direct towards an existing directory",
            "or file. Use '--reuse' to load previous results.",
            f"Provided path '{output_path}'",
        ]
    )
    if isfile(output_path):
        raise AssertionError(assert_msg + " leads to a file")
    if isdir(output_path):
        raise AssertionError(assert_msg + " leads to a directory.")
