"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from .database import DataBase, Record
from .region import GenomicRegion
import itertools
import logging
import os
from typing import IO, Iterator, List, Union


class ChromosomeWalker:
    """docstring for Walker"""

    __chromosome: bytes
    __db: DataBase
    __IH: IO

    def __init__(self, db: Union[DataBase, str], chromosome: bytes):
        super(ChromosomeWalker, self).__init__()
        self.__db = DataBase(db) if isinstance(db, str) else db
        if chromosome not in self.__db._chromosomes.keys():
            raise AssertionError
        self.__IH = open(
            os.path.join(self.__db._root, f"{chromosome.decode()}.bin"), "rb"
        )
        self.__chromosome = chromosome

    @property
    def db(self) -> DataBase:
        return self.__db

    def __jump_to_bin(self, start_from_nt: int = 0) -> None:
        """Move buffer pointer to the first record of a bin.

        Keyword Arguments:
            start_from_nt {int} -- position in nucleotides (default: {0})
        """
        position_in_bytes = self.__db._chromosomes.get_chromosome(
            self.__chromosome
        ).index[start_from_nt]
        self.__IH.seek(position_in_bytes)

    def read_next_record(self) -> bytes:
        """Read the next record.

        The buffer pointer moves to the start of the following record.

        Returns:
            bytes -- record bytes
        """
        return self.__IH.read(self.__db.record_byte_size)

    # def read_previous_record(self) -> bytes:
    #     """Read the previous record.

    #     The buffer pointer does not move.

    #     Returns:
    #         bytes -- record bytes
    #     """
    #     self.rewind()
    #     return self.read_next_record()

    def read_next_record_and_rewind(self) -> bytes:
        """Reads the next record.

        The buffer pointer does not move.

        Returns:
            bytes -- record bytes
        """
        record = self.read_next_record()
        self.rewind()
        return record

    def next_record_exists(self) -> bool:
        """Tries reading the next record.

        Returns:
            bool -- whether next record exists
        """
        return len(self.read_next_record_and_rewind()) != 0

    def rewind(self) -> None:
        """Rewind of one record.

        The buffer pointer moves to the beginning of the previous record.
        """
        self.__IH.seek(max(self.__IH.tell() - self.__db.record_byte_size, 0))

    # def skip(self) -> None:
    #     """Skip one record.

    #     The buffer pointer moves to the beginning of the following record.
    #     """
    #     self.__IH.seek(self.__IH.tell() + self.__db.record_byte_size)

    def fastforward(self, start_from_nt: int) -> None:
        """Jump to the first record at a given position.

        First the buffer pointer is moved to the beginning of the bin containing the
        start_from_nt position. Then, records are read an parsed until their start is
        greater than the start_from_nt value. Finally, the buffer pointer is moved to
        the beginning of the last parsed record.

        Arguments:
            IH {IO} -- database input handle
            chromosome {bytes} -- chromosome label
            start_from_nt {int} -- position to fastforward to (in nucleotides)
        """
        if start_from_nt == 0:
            self.__IH.seek(0)
            return
        self.__jump_to_bin(start_from_nt)

        record_start = 0
        while record_start < start_from_nt:
            if self.next_record_exists():
                record_start = Record(self.read_next_record(), self.__db.dtype)["start"]
            else:
                logging.warning(
                    " ".join(
                        [
                            "the specified location",
                            f"({self.__chromosome.decode()}:{start_from_nt})",
                            "is outside the database.",
                        ]
                    )
                )
                return

    def buffer(self, start_from_nt: int = 0, end_at_nt: int = -1) -> Iterator[Record]:
        """Buffer a chromosome's records.

        Buffer records from a chromosome within the specified region.
        To buffer the whole records, specify a [0, -1] region.

        Keyword Arguments:
            start_from_nt {int} -- region starting position (default: {0})
            end_at_nt {int} -- region end position (default: {-1})

        Yields:
            Iterator[Record] -- parsed record
        """
        self.fastforward(start_from_nt)
        record = self.read_next_record()
        while len(record) != 0:
            parsed_record = Record(record, self.__db.dtype)
            if parsed_record["start"] > end_at_nt and end_at_nt > 0:
                break
            yield parsed_record
            record = self.read_next_record()

    def walk_single_region(self, region: GenomicRegion) -> Iterator[List[Record]]:
        if region.chromosome != self.__chromosome:
            raise AssertionError
        focus_start, focus_end = region.focus
        record_list = list(self.buffer(focus_start, focus_end))
        yield record_list
        while region.can_increase_focus():
            region.increase_focus()
            new_focus_start, new_focus_end = region.focus
            record_list = list(
                itertools.chain(
                    *[
                        list(self.buffer(new_focus_start, focus_start)),
                        record_list,
                        list(self.buffer(focus_end, new_focus_end)),
                    ]
                )
            )
            yield record_list
            focus_start, focus_end = region.focus

    def walk_multiple_regions(self, region_set_list: List[List[GenomicRegion]]) -> None:
        for region_set_idx in range(len(region_set_list)):
            region_set = region_set_list[region_set_idx]
            for region_idx in range(len(region_set)):
                region = region_set[region_idx]
                for record_set in self.walk_single_region(region):
                    logging.info(
                        " ".join(
                            [
                                f"Window {region_set_idx}.{region_idx}",
                                f"({region.focus}):",
                                f"read {len(record_set)} records",
                            ]
                        )
                    )
