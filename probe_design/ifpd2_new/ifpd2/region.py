"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from .chromosome import ChromosomeData
import logging
from typing import List, Tuple


class GenomicRegion:
    """Details on genomic region with a central focus region."""

    __chrom: bytes
    __chromStart: int
    __chromEnd: int
    __focusStart: int
    __focusEnd: int
    __focusSize: int
    __focusStep: int

    def __init__(
        self,
        chrom_region: Tuple[bytes, int, int],
        focus_style: Tuple[float, float] = (1, 1),
    ):
        super(GenomicRegion, self).__init__()
        self.__init_region(*chrom_region)
        self.__init_focus(*focus_style)

    def __init_region(self, chrom: bytes, chromStart: int, chromEnd: int) -> None:
        """Initialize region coordinates.

        Arguments:
            chrom {bytes} -- chromosome name
            chromStart {int} -- chromosome start position (usually 0)
            chromEnd {int} -- chromosome end position
        """
        if len(chrom) == 0:
            raise AssertionError("chromosome cannot be empty")
        if chromStart < 0:
            raise AssertionError(
                f"start should be greater than or equal to 0: {chromStart}"
            )
        if chromEnd <= chromStart:
            raise AssertionError(
                f"end should be greater than start: {chromStart}-{chromEnd}"
            )
        self.__chrom = chrom
        self.__chromStart = chromStart
        self.__chromEnd = chromEnd

    def __init_focus_size(self, focus_style: float) -> None:
        """Initialize focus region size.

        Set the size of the central focus region. When focus style is lower than or
        equal to 1, it is considered to be a fraction of the genomic region size. When
        it is greater than 1, instead, it is considered to be the size in nt.

        Arguments:
            focus_style {float} -- focus region size in nt
                                   or as a fraction of the genomic region
        """
        if focus_style <= 0:
            raise AssertionError
        if focus_style > 1:
            if focus_style > self.__chromEnd - self.__chromStart:
                raise AssertionError
            self.__focusSize = int(focus_style)
        else:
            self.__focusSize = int((self.__chromEnd - self.__chromStart) * focus_style)

    def __init_focus_step(self, step_style: float) -> None:
        """Initialize step for focus region growth.

        Set the step at which the focus region grows. When step style is lower than or
        equal to 1, it is considered to be a fraction of the focus region size. When
        it is greater than 1, instead, it is considered to be in nt.

        Arguments:
            step_style {float} -- focus region growth step in nt
                                  or as a fraction of the focus region
        """
        if step_style <= 0:
            raise AssertionError
        if step_style > 1:
            self.__focusStep = int(step_style)
        else:
            self.__focusStep = int(self.__focusSize * step_style)

    def __init_focus(self, focus_style: float, focus_step_style: float) -> None:
        """Initialize focus region.

        More details on focus and focus step style in the __init_focus_size and
        __init_focus_step functions.

        Arguments:
            focus_style {float} -- focus region size in nt
                                   or as a fraction of the genomic region
            focus_step_style {float} -- focus region growth step in nt
                                        or as a fraction of the focus region
        """
        self.__init_focus_size(focus_style)
        self.__init_focus_step(focus_step_style)
        self.__focusStart = int(
            (self.__chromStart + self.__chromEnd) / 2 - self.__focusSize / 2
        )
        self.__focusEnd = int(
            (self.__chromStart + self.__chromEnd) / 2 + self.__focusSize / 2
        )

    @property
    def chromosome(self) -> bytes:
        return self.__chrom

    @property
    def start(self) -> int:
        return self.__chromStart

    @property
    def end(self) -> int:
        return self.__chromEnd

    @property
    def region(self) -> Tuple[int, int]:
        return (self.__chromStart, self.__chromEnd)

    @property
    def focus(self) -> Tuple[int, int]:
        return (self.__focusStart, self.__focusEnd)

    # @property
    # def focus_step(self) -> int:
    #     return self.__focusStep

    def can_increase_focus(self) -> bool:
        """Check if the focus region can grow more

        Returns:
            bool -- whether growth is possible
        """
        return self.focus != self.region

    def increase_focus(self) -> None:
        """Grow the focus region."""
        if self.can_increase_focus():
            logging.warning("cannot increase the focus region any further")
        self.__focusStart = max(
            self.__chromStart, self.__focusStart - self.__focusStep // 2
        )
        self.__focusEnd = min(self.__chromEnd, self.__focusEnd + self.__focusStep // 2)


class GenomicRegionBuilder:
    """docstring for GenomicRegionBuilder"""

    __chromosome: bytes
    __chromosome_size_nt: int
    __focus_style: Tuple[float, float]

    def __init__(
        self,
        chromosome_data: ChromosomeData,
        focus_style: float = 1,
        focus_step_style: float = 1,
    ):
        super(GenomicRegionBuilder, self).__init__()
        self.__chromosome = chromosome_data.name
        self.__chromosome_size_nt = chromosome_data.size_nt
        if focus_style <= 0:
            raise AssertionError
        if focus_step_style <= 0:
            raise AssertionError
        self.__focus_style = (focus_style, focus_step_style)

    def __build_overlapping(self, size: int, step: int) -> List[List[GenomicRegion]]:
        """Build lists of regions.

        Windows are generated based on size and step. Step should be smaller than size.
        Overlapping regions are placed in separate lists.

        Arguments:
            size {int} -- region size in nt
            step {int} -- region step in nt

        Returns:
            List[List[GenomicRegion]] -- generated region lists
        """
        if step >= size:
            raise AssertionError
        region_set_list: List[List[GenomicRegion]] = []
        for range_start in range(0, size, step):
            genomic_region_set: List[GenomicRegion] = []
            for start in range(range_start, self.__chromosome_size_nt, step):
                end = start + size
                if end < self.__chromosome_size_nt:
                    genomic_region_set.append(
                        GenomicRegion(
                            (self.__chromosome, start, end), self.__focus_style
                        )
                    )
            region_set_list.append(genomic_region_set)
        return region_set_list

    def __build_non_overlapping(
        self, size: int, step: int
    ) -> List[List[GenomicRegion]]:
        """Build a list of non-overlapping regions.

        Arguments:
            size {int} -- region size in nt
            step {int} -- region step in nt

        Returns:
            List[List[GenomicRegion]] -- generated region lists
        """
        if step != size:
            raise AssertionError
        genomic_region_set: List[GenomicRegion] = []
        for start in range(0, self.__chromosome_size_nt, step):
            end = start + size
            if end <= self.__chromosome_size_nt:
                genomic_region_set.append(
                    GenomicRegion((self.__chromosome, start, end), self.__focus_style)
                )
        return [genomic_region_set]

    def build_by_number(self, n: int) -> List[List[GenomicRegion]]:
        """Build a number of regions.

        Build a number of region instances covering the current chromosome.

        Arguments:
            n {int} -- number of regions

        Returns:
            List[List[GenomicRegion]] -- generated region lists
        """
        step: int = self.__chromosome_size_nt // n
        return self.build_by_size(step, step)

    def build_by_size(self, size: int, step_style: float) -> List[List[GenomicRegion]]:
        """Build regions by size and step (allows for overlap).

        Build region instances based on their size and step. Step is the distance
        between the start of a region and the start of the next region. Overlapping
        regions are placed in separate lists.

        Arguments:
            size {int} -- region size in nt
            step_style {float} -- region step in nt or as a fraction of the region size

        Returns:
            List[List[GenomicRegion]] -- generated region lists
        """
        step = int(step_style) if step_style > 1 else int(size * step_style)
        if step < size:
            return self.__build_overlapping(size, step)
        return self.__build_non_overlapping(size, step)
