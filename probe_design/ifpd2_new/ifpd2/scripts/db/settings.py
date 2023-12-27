"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from copy import copy
from ...const import DEFAULT_DATABASE_INDEX_BIN_SIZE
import logging
from os.path import isdir, isfile
from typing import Any, Dict, List, Set


class DBMakeSettings:
    _output_path: str
    _off_target_paths: Set[str]
    _melting_temperature_paths: Set[str]
    _secondary_structure_paths: Set[str]
    _bin_size: int = DEFAULT_DATABASE_INDEX_BIN_SIZE
    prefix: str = ""

    def __init__(self, output_path: str):
        super(DBMakeSettings, self).__init__()
        self.output_path = output_path
        self._off_target_paths = set()
        self._melting_temperature_paths = set()
        self._secondary_structure_paths = set()

    @property
    def output_path(self) -> str:
        return self._output_path

    @output_path.setter
    def output_path(self, output_path: str) -> None:
        if isdir(output_path) or isfile(output_path):
            raise AssertionError(f"'{output_path}' already exists")
        self._output_path = output_path

    @property
    def bin_size(self) -> int:
        return self._bin_size

    @bin_size.setter
    def bin_size(self, bin_size: int) -> None:
        if bin_size < 1:
            raise AssertionError
        self._bin_size = bin_size

    @property
    def off_target_paths(self) -> Set[str]:
        return copy(self._off_target_paths)

    @property
    def melting_temperature_paths(self) -> Set[str]:
        return copy(self._melting_temperature_paths)

    @property
    def secondary_structure_paths(self) -> Set[str]:
        return copy(self._secondary_structure_paths)

    def add_off_target_path(self, path: str) -> None:
        if not isfile(path):
            raise AssertionError(f"'{path}' not found.")
        self._off_target_paths.add(path)

    def add_off_target_path_list(self, path_list: List[str]) -> None:
        for path in path_list:
            self.add_off_target_path(path)

    def add_melting_temperature_path(self, path: str) -> None:
        if not isfile(path):
            raise AssertionError(f"'{path}' not found.")
        self._melting_temperature_paths.add(path)

    def add_melting_temperature_path_list(self, path_list: List[str]) -> None:
        for path in path_list:
            self.add_melting_temperature_path(path)

    def add_secondary_structure_path(self, path: str) -> None:
        if not isfile(path):
            raise AssertionError(f"'{path}' not found.")
        self._secondary_structure_paths.add(path)

    def add_secondary_structure_path_list(self, path_list: List[str]) -> None:
        for path in path_list:
            self.add_secondary_structure_path(path)

    def init_login(self) -> None:
        logging.info(f"output folder: '{self.output_path}'")
        logging.info(f"chromosome prefix: '{self.prefix}'")
        logging.info(f"hush: {self.off_target_paths}")
        logging.info(f"oligo-melting: {self.melting_temperature_paths}")
        logging.info(f"OligoArrayAux: {self.secondary_structure_paths}")

    def asdict(self) -> Dict[str, Any]:
        return dict(
            output_path=self.output_path,
            off_target_paths=self.off_target_paths,
            melting_temperature_paths=self.melting_temperature_paths,
            secondary_structure_paths=self.secondary_structure_paths,
            bin_size=self.bin_size,
            prefix=self.prefix,
        )
