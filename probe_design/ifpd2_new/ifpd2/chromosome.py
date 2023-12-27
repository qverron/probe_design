"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import copy
from . import const
from .io import get_dtype_length
import numpy as np  # type: ignore
import pandas as pd  # type: ignore
from rich.progress import Progress, TaskID  # type: ignore
from typing import Dict, List, Set, Tuple


class ChromosomeIndex:
    """ChromosomeIndex"""

    _bin_size: int
    _index: Dict[int, Tuple[float, float]]

    def __init__(self, bin_size: int):
        super(ChromosomeIndex, self).__init__()
        if bin_size < 1:
            raise AssertionError
        self._bin_size = bin_size

    def __init_index(self, chrom_db: pd.DataFrame) -> None:
        """Initialize chromosome index

        Creates an chromosome index with empty bins.

        Arguments:
            chrom_db {pd.DataFrame} -- chromosome database
        """
        self._index = {}
        chrom_size_nt = chrom_db["end"].values.max()
        for bin_id in range((chrom_size_nt // self._bin_size) + 1):
            self._index[bin_id] = (np.inf, 0)

    @property
    def bin_size(self) -> int:
        return self._bin_size

    @property
    def data(self) -> Dict[int, Tuple[float, float]]:
        return copy.copy(self._index)

    def __populate_bins(
        self,
        chrom_db: pd.DataFrame,
        record_byte_size: int,
        track: Tuple[Progress, TaskID],
    ) -> None:
        """Populate bins.

        Add extremities (byte positions) of each bin.

        Arguments:
            chrom_db {pd.DataFrame} -- chromosome database
            record_byte_size {int} -- size of a record in bytes (based on column dtype)
            track {Tuple[Progress, TaskID]} -- progress bar details
        """
        current_position = -1
        for i in range(chrom_db.shape[0]):
            track[0].update(track[1], advance=1)
            position_in_nt = chrom_db["start"].values[i]
            if position_in_nt <= current_position:
                raise AssertionError
            current_position = position_in_nt

            position_in_bytes = record_byte_size * i
            binned_to = position_in_nt // self._bin_size

            bin_start, bin_end = list(self._index[binned_to])
            bin_start = min(bin_start, position_in_bytes)
            bin_end = max(bin_end, position_in_bytes)
            self._index[binned_to] = (bin_start, bin_end)

    def __fill_empty_bins(self) -> None:
        """Populate empty bins of the index."""
        if not np.isfinite(self._index[0][0]):
            self._index[0] = (0, 0)
        for bin_id, (start, end) in self._index.items():
            if not np.isfinite(start):
                self._index[bin_id] = (
                    self._index[bin_id - 1][1],
                    self._index[bin_id - 1][1],
                )

    def build(
        self,
        chrom_db: pd.DataFrame,
        record_byte_size: int,
        track: Tuple[Progress, TaskID],
    ) -> None:
        """Build chromosome index.

        [description]

        Arguments:
            chrom_db {pd.DataFrame} -- chromosome database
            record_byte_size {int} -- size of a record in bytes (based on column dtype)
            track {Tuple[Progress, TaskID]} -- progress bar details
        """
        for colname in ("chromosome", "start", "end"):
            if colname not in chrom_db.columns:
                raise AssertionError(f"missing '{colname}' column")
        chromosome_set: Set[bytes] = set(chrom_db["chromosome"].values)
        if len(chromosome_set) != 1:
            raise AssertionError

        self.__init_index(chrom_db)
        self.__populate_bins(chrom_db, record_byte_size, track)
        self.__fill_empty_bins()

    def __getitem__(self, position_in_nt: int) -> int:
        if self._index is None:
            raise AssertionError
        binned_to = position_in_nt // self._bin_size
        if binned_to not in self._index:
            return -1

        position_in_bytes = int(self._index[binned_to][0])
        if not np.isfinite(position_in_bytes):
            return -1
        return position_in_bytes

    def __eq__(self, other) -> bool:
        condition = self.bin_size == other.bin_size
        condition = condition and self.data == other.data
        return condition


class ChromosomeData:
    """Contains information on a chromosome"""

    _name: bytes
    _size_nt: int
    _size_bytes: int
    _recordno: int
    _index: ChromosomeIndex
    _record_byte_size: int

    def __init__(
        self,
        chromosome_db: pd.DataFrame,
        dtype: Dict[str, str],
        index_bin_size: int,
        progress: Progress,
    ):
        super(ChromosomeData, self).__init__()

        if "chromosome" not in chromosome_db.columns:
            raise AssertionError
        selected_chrom = chromosome_db["chromosome"][0]
        if len(set(chromosome_db["chromosome"].values)) != 1:
            raise AssertionError

        self._record_byte_size = get_dtype_length(dtype)
        if self._record_byte_size <= 0:
            raise AssertionError

        self._name = selected_chrom
        self._recordno = chromosome_db.shape[0]
        self._size_nt = chromosome_db["end"].values.max()
        self._size_bytes = chromosome_db.shape[0] * self._record_byte_size

        self._build_index(chromosome_db, index_bin_size, progress)

    @property
    def name(self) -> bytes:
        return self._name

    @property
    def record_byte_size(self) -> int:
        return self._record_byte_size

    @property
    def size_nt(self) -> int:
        return self._size_nt

    @property
    def size_bytes(self) -> int:
        return self._size_bytes

    @property
    def recordno(self) -> int:
        return self._recordno

    @property
    def index(self) -> ChromosomeIndex:
        return copy.copy(self._index)

    def _build_index(
        self, chromosome_db: pd.DataFrame, index_bin_size: int, progress: Progress
    ) -> None:
        """Build a chromosome index.

        Build bin index based on chromosome length and bin size.

        Arguments:
            chromosome_db {pd.DataFrame} -- chromosome database
            index_bin_size {int} -- index bin size
            progress {Progress} -- progress instance for progress bar with rich
        """
        if index_bin_size <= 0:
            raise AssertionError
        indexing_track = progress.add_task(
            f"indexing {self._name.decode()}.bin",
            total=chromosome_db.shape[0],
            transient=True,
        )
        self._index = ChromosomeIndex(index_bin_size)
        self._index.build(
            chromosome_db, self._record_byte_size, (progress, indexing_track)
        )
        progress.remove_task(indexing_track)

    def __eq__(self, other) -> bool:
        condition = self.record_byte_size == other.record_byte_size
        condition = condition and self.size_nt == other.size_nt
        condition = condition and self.size_bytes == other.size_bytes
        condition = condition and self.recordno == other.recordno
        condition = condition and self.index == other.index
        return condition


class ChromosomeDict:
    """Wraps ChromosomeData instances into a dictionary
    and allow easier access to their attributes.
    """

    _index_bin_size: int
    _data: Dict[bytes, ChromosomeData]

    def __init__(self, index_bin_size: int = const.DEFAULT_DATABASE_INDEX_BIN_SIZE):
        super(ChromosomeDict, self).__init__()
        self._data = {}
        if index_bin_size <= 0:
            raise AssertionError
        self._index_bin_size = index_bin_size

    def __len__(self) -> int:
        return len(self._data)

    def keys(self) -> List[bytes]:
        return list(self._data.keys())

    def items(self) -> List[Tuple[bytes, ChromosomeData]]:
        return list(self._data.items())

    @property
    def index_bin_size(self) -> int:
        return self._index_bin_size

    @property
    def sizes_nt(self) -> Dict[bytes, int]:
        return {name: data.size_nt for name, data in self._data.items()}

    @property
    def sizes_bytes(self) -> Dict[bytes, int]:
        return {name: data.size_bytes for name, data in self._data.items()}

    @property
    def recordnos(self) -> Dict[bytes, int]:
        return {name: data.recordno for name, data in self._data.items()}

    def get_chromosome(self, chromosome: bytes) -> ChromosomeData:
        return copy.copy(self._data[chromosome])

    def add_chromosome(
        self, chromosome_db: pd.DataFrame, dtype: Dict[str, str], progress: Progress
    ) -> None:
        self._data[chromosome_db["chromosome"][0]] = ChromosomeData(
            chromosome_db, dtype, self._index_bin_size, progress
        )

    def __eq__(self, other) -> bool:
        condition = len(self) == len(other)
        condition = condition and self.items() == other.items()
        condition = condition and self.index_bin_size == other.index_bin_size
        return condition
