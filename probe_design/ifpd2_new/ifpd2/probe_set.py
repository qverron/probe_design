"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import itertools
import logging
import numpy as np  # type: ignore
import os
import pandas as pd  # type: ignore
from rich.progress import track  # type: ignore
import shutil
from typing import List

from .probe import OligoProbe


class OligoProbeSet:
    """Set of non-overlapping probes built from windows."""

    def __init__(self, probe_list):
        super(OligoProbeSet, self).__init__()
        self.probe_list = probe_list

    @property
    def probe_list(self):
        return self._probe_list

    @probe_list.setter
    def probe_list(self, probe_list):
        if not all(isinstance(p, OligoProbe) for p in probe_list):
            raise AssertionError
        self._probe_list = sorted(probe_list, key=lambda p: p.range[0])
        self._probe_tm_ranges = [p.tm_range for p in self.probe_list]
        self._tm_range = (
            np.min([t[0] for t in self._probe_tm_ranges]),
            np.max([t[1] for t in self._probe_tm_ranges]),
        )
        self._sizes = [p.size for p in self.probe_list]
        self._spreads = [p.spread for p in self.probe_list]
        self._probe_ranges = [p.range for p in self.probe_list]
        probe_starts = np.array([r[0] for r in self.probe_ranges])
        probe_ends = np.array([r[1] for r in self.probe_ranges])
        self._range = (probe_starts.min(), probe_ends.max())
        self._ds = 0 if len(probe_list) == 1 else probe_starts[1:] - probe_ends[:-1]
        self._d_mean = np.mean(self.ds)
        self._d_range = (np.min(self.ds), np.max(self.ds))
        tm = [p.data["Tm"].tolist() for p in self.probe_list]
        self._oligo_tm = list(itertools.chain(*tm))

    @property
    def range(self):
        return self._range

    @property
    def probe_ranges(self):
        return self._probe_ranges

    @property
    def ds(self):
        return self._ds

    @property
    def d_mean(self):
        return self._d_mean

    @property
    def d_range(self):
        return self._d_range

    @property
    def tm_range(self):
        return self._tm_range

    @property
    def probe_tm_ranges(self):
        return self._probe_tm_ranges

    @property
    def oligo_tm(self):
        return self._oligo_tm

    @property
    def sizes(self):
        return self._sizes

    @property
    def spreads(self):
        return self._spreads

    @property
    def featDF(self):
        mean_score = (1 / 3) * sum(
            [
                np.std(self.sizes) / np.mean(self.sizes),
                np.std(self.spreads) / np.mean(self.spreads),
                np.std(self.oligo_tm) / np.mean(self.oligo_tm),
            ]
        )
        df = pd.DataFrame(
            [
                self.range[0],
                self.range[1],
                len(self.probe_list),
                np.diff(self.range)[0],
                np.mean(self.sizes),
                np.std(self.sizes),
                np.min(self.spreads),
                np.max(self.spreads),
                np.mean(self.spreads),
                np.std(self.spreads),
                self.d_range[0],
                self.d_range[1],
                self.d_mean,
                np.std(self.ds),
                np.diff(self.tm_range)[0],
                np.mean(self.oligo_tm),
                np.std(self.oligo_tm),
                mean_score,
            ]
        ).transpose()
        df.columns = [
            "start",
            "end",
            "nProbes",
            "size",
            "size_mean",
            "size_std",
            "spread_min",
            "spread_max",
            "spread_mean",
            "spread_std",
            "d_min",
            "d_max",
            "d_mean",
            "d_std",
            "tm_range",
            "tm_mean",
            "tm_std",
            "score",
        ]
        return df

    def export(self, path):
        if os.path.isfile(path):
            raise AssertionError
        if os.path.isdir(path):
            shutil.rmtree(path)
        os.mkdir(path)

        self.featDF.to_csv(os.path.join(path, "set.tsv"), "\t", index=False)

        pd.concat([p.featDF for p in self.probe_list]).reset_index(drop=True).to_csv(
            os.path.join(path, "probes.tsv"), "\t"
        )

        with open(os.path.join(path, "set.bed"), "w+") as BH:
            for pi in range(len(self.probe_list)):
                probe = self.probe_list[pi]
                probe.export(os.path.join(path, f"probe_{pi}"))
                for i in range(probe.data.shape[0]):
                    oligo = probe.data.iloc[i]
                    BH.write(
                        "".join(
                            [
                                f">probe_{pi}:{oligo['name']}",
                                f":{oligo['chromosome']}:",
                                f"{oligo['start']}-{oligo['end']}\n",
                            ]
                        )
                    )
                    BH.write(f"{probe.data.iloc[i]['sequence']}\n")


class OligoProbeSetBuilder:
    """Class to build OligoProbeSet objects."""

    probe_set_list: List = []

    def __init__(self, out_path):
        super(OligoProbeSetBuilder, self).__init__()
        if os.path.isfile(out_path):
            raise AssertionError
        self.out_path = out_path
        if os.path.isdir(out_path):
            shutil.rmtree(out_path)
        os.mkdir(out_path)

    @staticmethod
    def __build_probe_set_list(window_list, i):
        probe_set_list = [(p,) for p in window_list[i]]

        for w in window_list[(i + 1) :]:
            if len(w) == 0:
                continue

            current_probe_set_list = list(tuple(probe_set_list))
            probe_set_list = []

            for probe in w:
                for probe_set in current_probe_set_list:
                    current_probe_set = list(probe_set)
                    current_probe_set.append(probe)
                    probe_set_list.append(current_probe_set)

        return probe_set_list

    def build(self, probe_candidates):
        for (wSet, window_list) in probe_candidates.items():
            window_list = list(window_list.values())

            non_empty_ids = [
                i for i in range(len(window_list)) if len(window_list[i]) != 0
            ]

            if not non_empty_ids:
                logging.warning(f"No probe candidates found, dropped. [ws{wSet+1}]")
                continue
            i = non_empty_ids[0]

            probe_set_list = self.__build_probe_set_list(window_list, i)

            self.probe_set_list.extend(probe_set_list)
            logging.info(
                "".join(
                    [
                        f"Built {len(probe_set_list)} probe set candidates ",
                        f"from window set #{wSet+1}",
                    ]
                )
            )

        self.probe_set_list = [OligoProbeSet(ps) for ps in self.probe_set_list]
        logging.info(
            f"Built {len(self.probe_set_list)} " + "probe set candidates in total."
        )

        if self.probe_set_list:
            pd.concat([ps.featDF for ps in self.probe_set_list]).reset_index(
                drop=True
            ).to_csv(os.path.join(self.out_path, "probe_sets.tsv"), "\t")

    def export(self):
        logging.info("Exporting probe sets...")
        for psi in track(range(len(self.probe_set_list)), description="Probe set"):
            probe_set = self.probe_set_list[psi]
            probe_set.export(os.path.join(self.out_path, f"probe_set_{psi}"))
