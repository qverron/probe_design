"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import click  # type: ignore
import copy
from ...const import CONTEXT_SETTINGS, DEFAULT_DATABASE_INDEX_BIN_SIZE
from ...const import dtype_header_features, dtype_sequence_features
from ...const import database_columns, dtype_hush, dtype_melting, dtype_secondary
from ... import asserts as ass
from ...chromosome import ChromosomeDict
from ...io import parse_hush, parse_melting, parse_secondary
from .settings import DBMakeSettings
import logging
import numpy as np  # type: ignore
import os
import pandas as pd  # type: ignore
import pickle
from rich.progress import Progress, track  # type: ignore
from typing import Callable, Dict, List, Set, Tuple


@click.command(
    name="make",
    context_settings=CONTEXT_SETTINGS,
    help="""
Assembles files into a database. At least one of the following is required,
to retrieve sequence: hush output (-O), or oligo-melting output (-T).

The name of each record should be in the format: 'name pos=chrom:start-end'.

Accepted formats:
hush            single-line sequence fasta file, where a ", XXX" is appended after the
                  sequence. XXX is the number of detected off-targets.
oligo-melting   .tsv file with six (6) columns: name, dG, dH, dS, Tm, Seq.
                  NOTE: the first line of the file is skipped during parsing because
                  it is expected to be a header line.
OligoArrayAux   .ct (Connectivity Table) file. Details on format available here:
                  https://rna.urmc.rochester.edu/Text/File_Formats.html

Files are merged by: hush header, oligo-melting name, and
OligoArrayAux name (3rd field of 1st column).
""",
)
@click.argument("output", nargs=1, type=click.Path(exists=False))
@click.option(
    "--hush",
    "-O",
    metavar="hush_path",
    nargs=1,
    multiple=True,
    type=click.Path(exists=True),
    help="Path to hush output. Format details above.",
)
@click.option(
    "--melting",
    "-T",
    metavar="hush_oligomelting_pathpath",
    nargs=1,
    multiple=True,
    type=click.Path(exists=True),
    help="Path to oligo-melting output. Format details above.",
)
@click.option(
    "--secondary",
    "-S",
    metavar="oligoarrayaux_path",
    nargs=1,
    multiple=True,
    type=click.Path(exists=True),
    help="Path to OligoArrayAux output. Format details above.",
)
@click.option(
    "--prefix",
    "-p",
    type=click.STRING,
    default="",
    help="Prefix to be added to chromosome labels.",
)
@click.option(
    "--binsize",
    "-b",
    default=DEFAULT_DATABASE_INDEX_BIN_SIZE,
    type=click.INT,
    help=f"Index bin size. Default: {DEFAULT_DATABASE_INDEX_BIN_SIZE}",
)
def main(
    output: str,
    hush: List[str],
    melting: List[str],
    secondary: List[str],
    prefix: str,
    binsize: int,
) -> None:
    settings = DBMakeSettings(output)
    settings.add_off_target_path_list(hush)
    settings.add_melting_temperature_path_list(melting)
    settings.add_secondary_structure_path_list(secondary)
    settings.bin_size = binsize
    settings.prefix = prefix

    if (
        len(settings.off_target_paths) == 0
        and len(settings.melting_temperature_paths) != 0
    ):
        raise AssertionError("please provide either --hush or --melting")

    os.mkdir(settings.output_path)
    dbdf = pd.DataFrame(columns=["name"])
    dbdf.set_index("name", inplace=True)

    dbdf = populate_db(dbdf, settings.off_target_paths, parse_hush, "hush")
    dbdf = populate_db(
        dbdf, settings.melting_temperature_paths, parse_melting, "oligo-melting"
    )
    dbdf = populate_db(
        dbdf, settings.secondary_structure_paths, parse_secondary, "OligoArrayAux"
    )

    dbdf = reduce_sequence_columns(dbdf)
    dbdf, dtype_sequence = parse_sequences(dbdf)
    dbdf, dtype_header = parse_record_headers(dbdf, prefix)

    dtype = {}
    dtype.update(dtype_melting)
    dtype.update(dtype_hush)
    dtype.update(dtype_secondary)
    dtype.update(dtype_sequence)
    dtype.update(dtype_header)

    for column in dtype:
        if column not in dbdf.columns:
            dbdf["column"] = np.repeat(np.nan, dbdf.shape[0])
    dbdf = dbdf.loc[:, database_columns]

    write_database(dbdf, dtype, settings)

    logging.info("Done. :thumbs_up: :smiley:")


def parse_input(
    path_list: Set[str], parse_function: Callable, software_name: str
) -> pd.DataFrame:
    if not path_list:
        return pd.DataFrame()
    logging.info(f"parsing {software_name} output")
    return pd.concat([parse_function(path) for path in path_list])


def populate_db(
    db: pd.DataFrame,
    path_list: Set[str],
    parse_function: Callable,
    software_name: str,
) -> pd.DataFrame:
    if not path_list:
        return db
    parsed_db = parse_input(path_list, parse_function, software_name)
    logging.info(f"adding {software_name} output to database")
    return db.merge(
        parsed_db,
        how="outer",
        left_index=True,
        right_index=True,
    )


def reduce_sequence_columns(df: pd.DataFrame) -> pd.DataFrame:
    logging.info("discarding redundant sequence columns")
    seq_columns = [c for c in df.columns if "sequence" in c]
    if len(seq_columns) == 1:
        return df
    df.drop(seq_columns[1:], axis=1, inplace=True)
    df.rename(columns={seq_columns[0]: "sequence"}, inplace=True)
    return df


def parse_sequences(df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, str]]:
    logging.info("adding sequence feature columns: length and GC-content")
    sequence_length_list: List[int] = []
    gc_content_list: List[float] = []

    for sequence in track(
        df["sequence"].values, description="calculating GC-content", transient=True
    ):
        sequence = sequence.upper()
        sequence_length_list.append(len(sequence))
        gc_content_list.append(
            (sequence.count(b"G") + sequence.count(b"C")) / len(sequence)
        )

    df["gc_content"] = gc_content_list
    dtype = copy.copy(dtype_sequence_features)
    dtype["sequence"] = f"|S{max(sequence_length_list)}"
    return (df.astype(dtype), dtype)


def parse_record_headers(
    db: pd.DataFrame, chromosome_prefix: str = ""
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    logging.info("adding header feature columns: name, chromosome, start, end")
    name_list: List[str] = []
    name_length_set: Set[int] = set()
    chromosome_list: List[str] = []
    chromosome_length_set: Set[int] = set()
    start_list: List[int] = []
    end_list: List[int] = []

    for record in track(
        db.itertuples(),
        total=db.shape[0],
        description="parsing record headers",
        transient=True,
    ):
        name, position = record.Index.split(" ")
        name_list.append(name)
        name_length_set.add(len(name))
        chromosome, extremes = position.split("=")[1].split(":")
        chromosome_list.append(f"{chromosome_prefix}{chromosome}")
        chromosome_length_set.add(len(f"{chromosome_prefix}{chromosome}"))
        start, end = [int(x) for x in extremes.split("-")]
        if (end - start) != len(record.sequence):
            raise AssertionError(f"{end - start} != {len(record.sequence)}")
        start_list.append(start)
        end_list.append(end)

    db["name"] = name_list
    db["chromosome"] = chromosome_list
    db["start"] = start_list
    db["end"] = end_list
    ass.ert_in_dtype(db["start"].values.max(), "u4")
    ass.ert_in_dtype(db["end"].values.max(), "u4")
    db.reset_index(drop=True, inplace=True)

    dtype = copy.copy(dtype_header_features)
    dtype["name"] = f"|S{max(name_length_set)}"
    dtype["chromosome"] = f"|S{max(chromosome_length_set)}"
    return (db.astype(dtype), dtype)


def write_database(
    dbdf: pd.DataFrame, dtype: Dict[str, str], settings: DBMakeSettings
) -> None:
    with Progress() as progress:
        chromosome_set: Set[bytes] = set(dbdf["chromosome"].values)
        chromosome_data = ChromosomeDict(settings.bin_size)
        chromosome_task = progress.add_task(
            "exporting chromosome",
            total=len(chromosome_set),
            transient=True,
        )
        for selected_chrom in chromosome_set:
            chromosome_db = dbdf.loc[selected_chrom == dbdf["chromosome"], :]

            logging.info(f"sorting records for {selected_chrom.decode()}")
            chromosome_db.sort_values(
                by="start", axis=0, kind="mergesort", inplace=True
            )
            logging.info(f"building index for {selected_chrom.decode()}")
            chromosome_data.add_chromosome(chromosome_db, dtype, progress)

            with open(
                os.path.join(settings.output_path, f"{selected_chrom.decode()}.bin"),
                "wb",
            ) as IH:
                writing_track = progress.add_task(
                    f"writing {selected_chrom.decode()}.bin",
                    total=chromosome_db.shape[0],
                    transient=True,
                )
                for record in chromosome_db.to_records(
                    index=False, column_dtypes=dtype
                ):
                    IH.write(record.tobytes())
                    progress.update(writing_track, advance=1)
            progress.update(chromosome_task, advance=1)

    logging.info("writing db.pickle")
    with open(os.path.join(settings.output_path, "db.pickle"), "wb") as OH:
        pickle.dump(
            dict(chromosomes=chromosome_data, dtype=dtype, args=settings.asdict()), OH
        )
