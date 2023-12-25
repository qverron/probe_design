"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from dataclasses import dataclass
from os.path import isdir, isfile
from typing import Optional, Tuple


def path_exists(path: str) -> bool:
    return isdir(path) or isfile(path)


@dataclass(frozen=True)
class Folder:
    path: str
    exists: bool = False

    def __post_init__(self):
        if (
            self.exists
            and not isdir(self.path)
            or not self.exists
            and path_exists(self.path)
        ):
            raise AssertionError((self.path, self.exists))


@dataclass(frozen=True)
class File:
    path: str
    exists: bool = False

    def __post_init__(self):
        if (
            self.exists
            and not isfile(self.path)
            or not self.exists
            and path_exists(self.path)
        ):
            raise AssertionError


@dataclass(frozen=True)
class PositiveInteger:
    n: int

    def __post_init__(self):
        if self.n < 1:
            raise AssertionError


@dataclass(frozen=True)
class NonNegativeFloat:
    n: float

    def __post_init__(self):
        if self.n < 0:
            raise AssertionError


@dataclass(frozen=True)
class NonNegativeInteger:
    n: int

    def __post_init__(self):
        if self.n <= 0:
            raise AssertionError


@dataclass(frozen=True)
class PositiveFloat:
    n: float
    limit: Optional[float] = None
    limit_included: bool = True

    def __post_init__(self):
        if self.limit is None:
            if self.n <= 0:
                raise AssertionError

        elif self.limit_included:
            if not 0 < self.n <= self.limit:
                raise AssertionError
        elif not 0 < self.n < self.limit:
            raise AssertionError


@dataclass(frozen=True)
class GenomicRegion:
    start: int
    end: int

    def __post_init__(self):
        if self.start < 0:
            raise AssertionError
        if self.end < self.start and self.end != -1:
            raise AssertionError

    def astuple(self) -> Tuple[int, int]:
        return (self.start, self.end)


@dataclass(frozen=True)
class NonNegativeIntInterval:
    start: int
    end: int

    def __post_init__(self):
        if self.start < 0:
            raise AssertionError
        if self.end < self.start:
            raise AssertionError

    def astuple(self) -> Tuple[int, int]:
        return (self.start, self.end)


@dataclass(frozen=True)
class QueryWindow:
    size: Optional[int]
    shift: Optional[float]

    def __post_init__(self):
        if self.size is not None and self.size < 1:
            raise AssertionError
        if not 0 < self.shift <= 1:
            raise AssertionError

    def astuple(self) -> Tuple[Optional[int], Optional[float]]:
        return (self.size, self.shift)


@dataclass(frozen=True)
class QueryFocus:
    size: float
    step: float

    def __post_init__(self):
        if self.size <= 0:
            raise AssertionError
        if self.step <= 0:
            raise AssertionError


@dataclass(frozen=True)
class FreeEnergyInterval:
    start: float
    end: float

    def astuple(self) -> Tuple[float, float]:
        return (self.start, self.end)


@dataclass(frozen=True)
class GCRange:
    low: float
    high: float

    def __post_init__(self):
        if not 0 <= self.low <= 1:
            raise AssertionError
        if not 0 <= self.high <= 1:
            raise AssertionError
        if self.low > self.high:
            raise AssertionError
