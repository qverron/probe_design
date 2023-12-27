"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from . import asserts as ass
from .dataclasses import GCRange
from . import loggingg
import numpy as np  # type: ignore
import logging
import pandas as pd  # type: ignore
from tqdm import tqdm  # type: ignore
from typing import List, Tuple


class Oligo:
    """Oligo database line/record parser.
    Presents oligo values as properties.
    """

    _data = None

    __colnames = [
        "name",
        "chrom",
        "start",
        "end",
        "tm_dG",
        "dH",
        "dS",
        "Tm",
        "seq",
        "nOT",
        "k",
        "ss_dG",
        "GC",
    ]

    def __init__(self, oligo, i):
        super(Oligo, self).__init__()
        ass.ert_type(i, int, "oligo id")
        if i < 0:
            raise AssertionError
        self._raw_data = oligo.strip().split("\t")
        for i in [2, 3, 9, 10]:
            self._raw_data[i] = int(self._raw_data[i])
        for i in [4, 5, 6, 7, 11, 12]:
            self._raw_data[i] = float(self._raw_data[i])
        self._idx = i

    @property
    def idx(self):
        return self._idx

    @property
    def raw_data(self):
        return self._raw_data

    @property
    def data(self):
        if isinstance(self._data, type(None)):
            self._data = pd.DataFrame(
                self._raw_data, columns=[self._idx], index=self.__colnames
            ).transpose()
        return self._data

    @property
    def start(self):
        return self.data["start"].values[0]

    @property
    def end(self):
        return self.data["end"].values[0]

    @property
    def score(self):
        if "score" in self.data.columns:
            return self.data["score"].values[0]
        return np.nan

    @staticmethod
    def __norm_value_in_range(v, r):
        if r[1] == 0:
            return np.nan
        return (v - r[0]) / (r[1] - r[0])

    def __update_score_by_nOT(self, F):
        nOT = self.data["nOT"].values
        if nOT <= F[0]:
            self.data["score"] = 0
            return None
        if nOT > F[1]:
            self.data["score"] = np.inf
            return None
        return self.__norm_value_in_range(nOT, F)

    def __update_score_by_dG_Tm(self, Gs):
        ss_dG = self.data["ss_dG"].values
        tm_dG = self.data["tm_dG"].values
        if ss_dG >= tm_dG * min(Gs):
            self.data["score"] = 0
            return None
        if ss_dG < tm_dG * max(Gs):
            self.data["score"] = np.inf
            return None
        return self.__norm_value_in_range(ss_dG, [tm_dG * f for f in Gs])

    def __update_score_by_dG_Gs(self, Gs):
        ss_dG = self.data["ss_dG"].values
        if ss_dG >= Gs[0]:
            self.data["score"] = 0
            return None
        if ss_dG < Gs[1]:
            self.data["score"] = np.inf
            return None
        return self.__norm_value_in_range(ss_dG, Gs)

    def add_score(self, F, Gs):
        norm_nOT = self.__update_score_by_nOT(F)
        if norm_nOT is None:
            return
        if all(x >= 0 for x in Gs):
            norm_ss_dG = self.__update_score_by_dG_Tm(Gs)
        else:
            norm_ss_dG = self.__update_score_by_dG_Gs(Gs)
        if norm_ss_dG is None:
            return
        self.data["score"] = np.mean([norm_nOT, norm_ss_dG])


class OligoGroup:
    """Allows to select oligos from a group based on a "focus" window of
    interest. The window can be expanded to the closest oligo or to retain at
    least a given number of oligos.
    """

    _focus_window = None  # Left-close, right-open
    _oligos_in_focus_window = None
    _oligos_passing_score_filter = None

    def __init__(self, oligos, logger=logging.getLogger()):
        super(OligoGroup, self).__init__()
        self.logger = logger
        self._data = pd.concat(
            [o.to_dataframe() for o in oligos if np.isfinite(o["score"])],
            ignore_index=True,
        )
        self._data = self._data.loc[self._data["score"] <= 1, :]
        self._oligos_passing_score_filter = self._data.shape[0]

    @property
    def data(self):
        return self._data.copy()

    @property
    def focus_window(self):
        return self._focus_window

    @focus_window.setter
    def focus_window(self, focus_window):
        ass.ert_type(focus_window, tuple, "focus window")
        if len(focus_window) != 2:
            raise AssertionError
        if focus_window[1] <= focus_window[0]:
            raise AssertionError
        self._focus_window = focus_window

    @property
    def focus_window_repr(self):
        return f"[{self.focus_window[0]}:{self.focus_window[1]})"

    @property
    def focus_window_size(self):
        return self.focus_window[1] - self.focus_window[0]

    @property
    def oligos_in_focus_window(self):
        return self._oligos_in_focus_window

    @property
    def oligos_passing_score_filter(self):
        return self._oligos_passing_score_filter

    @property
    def usable_oligos(self):
        if isinstance(self.oligos_in_focus_window, type(None)):
            return self.oligos_passing_score_filter
        return np.logical_and(
            self.oligos_in_focus_window, self.oligos_passing_score_filter
        )

    def get_focused_oligos(self, onlyUsable=False):
        if onlyUsable:
            return self._data.loc[self.usable_oligos]
        return self._data.loc[self.oligos_in_focus_window, :]

    def get_n_focused_oligos(self, onlyUsable=False):
        if self.usable_oligos.sum() == 0:
            return 0
        return self.get_focused_oligos(onlyUsable).shape[0]

    def focus_all(self):
        # Use all oligos to define a focus region
        self.set_focus_window(
            self._data.loc[:, "start"].min(), self._data.loc[:, "start"].max() + 1
        )

    def set_focus_window(self, start, end, verbose=True):
        # Set a sub-window of interest to focus on
        self.focus_window = (start, end)

        start_condition = self._data["start"].values >= self.focus_window[0]
        end_condition = self._data["end"].values < self.focus_window[1]
        self._oligos_in_focus_window = np.logical_and(start_condition, end_condition)

        if verbose:
            nOligos = self.get_n_focused_oligos()
            nOligosUsable = self.get_n_focused_oligos(True)
            self.logger.info(
                " ".join(
                    [
                        f"Set focus region to {self.focus_window_repr}",
                        f"({nOligos} oligos, {nOligosUsable} usable)",
                    ]
                )
            )

    def expand_focus_to_n_oligos(self, n, verbose=True):
        # Expand the sub-window of interest to retrieve at least n oligos
        if isinstance(self.focus_window, type(None)):
            raise AssertionError

        if n <= self.get_n_focused_oligos():
            return

        for _ in range(n - self.get_n_focused_oligos()):
            if not self.expand_focus_to_closest():
                return

        if verbose:
            nOligos = self.get_n_focused_oligos()
            nOligosUsable = self.get_n_focused_oligos(True)
            self.logger.info(
                " ".join(
                    [
                        f"Expanded focus region to {self.focus_window}",
                        f" ({nOligos} oligos, {nOligosUsable} usable)",
                    ]
                )
            )

    def expand_focus_by_step(self, step, verbose=True):
        # Expand the current focus window of a given step (in nt)
        if step >= 0:
            raise AssertionError

        if (
            self.focus_window[0] <= self._data["start"].min()
            and self.focus_window[1] >= self._data["end"].max()
        ):
            self.logger.warning(
                " ".join(
                    [
                        "Cannot expand the focus region any further ",
                        "(all oligos already included)",
                    ]
                )
            )
            return False

        new_focus_start, new_focus_end = self.focus_window
        new_focus_start -= step / 2
        new_focus_end += step / 2
        new_focus_start = np.max([new_focus_start, self._data["start"].min()])
        new_focus_end = np.min([new_focus_end, self._data["end"].max()])

        self.set_focus_window(new_focus_start, new_focus_end, verbose)
        return True

    def __get_focus_extremes(self):
        earl = self._data["start"].values < self.focus_window[0]
        max_start = self._data.loc[earl, "start"].max()
        d_earl = self.focus_window[0] - max_start if earl.sum() != 0 else np.inf
        late = self._data["end"].values >= self.focus_window[1]
        min_end = self._data.loc[late, "end"].min()
        d_late = min_end - self.focus_window[1] if late.sum() != 0 else np.inf
        return (max_start, min_end, d_earl, d_late)

    def __set_focus_windows_for_inifinite_extreme(self):
        max_start, min_end, d_earl, d_late = self.__get_focus_extremes()
        if np.isinf(d_late):
            self.set_focus_window(max_start, self.focus_window[1], False)
        if np.isinf(d_earl):
            self.set_focus_window(self.focus_window[0], min_end + 1, False)

    def expand_focus_to_closest(self):
        # Expand the sub-window of interest to add the closest oligo
        # Return False if not possible (e.g., all oligos already included)
        if self.get_n_focused_oligos() == self._data.shape[0]:
            self.logger.warning(
                " ".join(
                    [
                        "Cannot expand the focus region any further",
                        "(all oligos already included)",
                    ]
                )
            )
            return False

        max_start, min_end, d_earl, d_late = self.__get_focus_extremes()

        if np.isinf(d_late) and np.isinf(d_earl):
            return False

        if np.isinf(d_late) or np.isinf(d_earl):
            self.__set_focus_windows_for_inifinite_extreme()
            return True

        if d_earl <= d_late:
            self.set_focus_window(max_start, self.focus_window[1], False)
        else:
            self.set_focus_window(self.focus_window[0], min_end + 1, False)

        return True

    def apply_threshold(self, threshold):
        # Unfocuses oligos with score higher than the threshold
        if threshold > 1 or threshold < 0:
            raise AssertionError
        self._oligos_passing_score_filter = self._data["score"] <= threshold

    def reset_threshold(self):
        self.apply_threshold(1)

    def discard_focused_oligos_safeDist(self, safeDist):
        # Discard focused oligos that are not within a safe distance from the
        # CFR borders
        start = self.data.loc[self.oligos_in_focus_window, "start"].min()
        end = self.data.loc[self.oligos_in_focus_window, "end"].max() + 1
        self.__discard_oligos_in_range(start + safeDist, end - safeDist)

    def __check_oligos_to_discard_safeN(self, safeN, start, end, D):
        oData = self.get_focused_oligos()
        c = 1
        while c <= safeN:
            passing_oData = oData.loc[oData["start"] > start, :]
            if passing_oData.shape[0] == 0:
                self.logger.info("Not enough oligos, skipped discard step.")
                return False
            start = passing_oData["start"].values[0]
            start += len(passing_oData["sequence"].values[0]) + D
            c += 1

        c = 1
        while c <= safeN:
            passing_oData = oData.loc[oData["end"] <= end]
            if passing_oData.shape[-1] == 0:
                self.logger.info("Not enough oligos, skipped discard step.")
                return False
            end = passing_oData["end"].values[-1]
            end -= len(passing_oData["sequence"].values[-1]) - D
            c += 1

        return True

    def discard_focused_oligos_safeN(self, safeN, D):
        # Discard focused oligos that are neither the first nor the last safeN
        start = self.data.loc[self.oligos_in_focus_window, "start"].values[0]
        start += (
            len(self.data.loc[self.oligos_in_focus_window, "sequence"].values[0]) + D
        )
        end = self.data.loc[self.oligos_in_focus_window, "end"].values[-1]
        end += (
            len(self.data.loc[self.oligos_in_focus_window, "sequence"].values[-1]) + D
        )

        if not self.__check_oligos_to_discard_safeN(safeN, start, end, D):
            return

        if not self.__discard_oligos_in_range(start, end):
            self.logger.info(f"No oligos to discard in range [{start}:{end}).")

    def __discard_oligos_in_range(self, start, end):
        if start >= end:
            return False
        start_condition = self.data["end"] < start
        end_condition = self.data["start"] >= end
        keep_condition = np.logical_or(start_condition, end_condition)
        nDiscarded = self.oligos_in_focus_window.sum() - keep_condition.sum()
        self.__apply_keep_condition(keep_condition)
        self.logger.info(
            f"Discarded {nDiscarded}" + f" oligos from the [{start}:{end}) range."
        )
        return True

    def __apply_keep_condition(self, keep_condition):
        self._data = self._data.loc[keep_condition, :]
        self._oligos_in_focus_window = self._oligos_in_focus_window[keep_condition]
        self._oligos_passing_score_filter = self._oligos_passing_score_filter[
            keep_condition
        ]


def select_by_GC(
    oligos_list: List[Tuple[str, str]],
    kmer_size: int,
    gc_range: GCRange = GCRange(0.35, 0.85),
) -> List[Tuple[str, str]]:
    return [
        (header, sequence)
        for header, sequence in tqdm(oligos_list, desc="GC check", leave=False)
        if (sequence.count("C") + sequence.count("G")) / kmer_size >= gc_range.low
        and (sequence.count("C") + sequence.count("G")) / kmer_size <= gc_range.high
    ]
