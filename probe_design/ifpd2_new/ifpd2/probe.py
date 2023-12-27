"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import configparser as cp
import numpy as np  # type: ignore
import os
import pandas as pd  # type: ignore
import shutil

from . import asserts as ass


class OligoProbe:
    """Converts a DataFrame of oligo data into an OligoProbe."""

    def __init__(self, oligo_data):
        super(OligoProbe, self).__init__()
        self.data = oligo_data

    @property
    def data(self):
        return self._data.copy()

    @data.setter
    def data(self, oligo_data):
        if not isinstance(oligo_data, pd.DataFrame):
            raise AssertionError
        required_columns = ["start", "end", "Tm"]
        if any(col not in oligo_data.columns for col in required_columns):
            raise AssertionError
        self._data = oligo_data
        self._range = (self._data["start"].min(), self._data["end"].max())
        self._size = self._range[1] - self._range[0]
        oDists = self._data["start"].values[1:] - self._data["end"].values[:-1]
        self._spread = np.std(oDists)
        self._d_range = (oDists.min(), oDists.max())
        self._d_mean = oDists.mean()
        self._tm_range = (self._data["Tm"].min(), self._data["Tm"].max())

    @property
    def n_oligos(self):
        return self._data.shape[0]

    @property
    def path(self):
        return self._data.index

    @property
    def range(self):
        return self._range

    @property
    def tm_range(self):
        return self._tm_range

    @property
    def d_range(self):
        return self._d_range

    @property
    def d_mean(self):
        return self._d_mean

    @property
    def spread(self):
        return self._spread

    @property
    def size(self):
        return self._size

    @property
    def featDF(self):
        df = pd.DataFrame(
            [
                self.range[0],
                self.range[1],
                self.n_oligos,
                self.size,
                self.spread,
                self.d_range[0],
                self.d_range[1],
                self.d_mean,
                np.diff(self.tm_range)[0],
            ]
        ).transpose()
        df.columns = [
            "start",
            "end",
            "nOligos",
            "size",
            "spread",
            "d_min",
            "d_max",
            "d_mean",
            "tm_range",
        ]
        return df

    def __repr__(self):
        rep = f"<OligoProbe[{self.range[0]}:{self.range[1]}"
        rep += f":{self.size}:{self.spread}]>"
        return rep

    def count_shared_oligos(self, probe):
        # Counts oligos shared with another probe
        # based on their paths
        return np.intersect1d(self.path, probe.path).shape[0]

    def export(self, path):
        if os.path.isfile(path):
            raise AssertionError
        if os.path.isdir(path):
            shutil.rmtree(path)
        os.mkdir(path)

        self.featDF.to_csv(os.path.join(path, "probe.tsv"), "\t", index=False)
        self.data.to_csv(os.path.join(path, "oligos.tsv"), "\t", index=False)

        with open(os.path.join(path, "probe.fa"), "w+") as BH:
            for i in range(self.data.shape[0]):
                oligo = self.data.iloc[i]
                BH.write(
                    "".join(
                        [
                            f">{oligo['name']}:",
                            f"{oligo['chromosome']}:{oligo['start']}-{oligo['end']}\n",
                        ]
                    )
                )
                BH.write(f"{oligo['sequence']}\n")


class OligoPathBuilder:
    """OligoPathBuilder contains methods to select non-overlapping sets of
    oligos (i.e., probes) from an oligo DataFrame in input.
    """

    N = int(48)  # Number of oligos per probe
    D: int = int(2)  # Minimum distance between consecutive oligos
    Tr = 10.0  # Melting temperature range half-width
    Ps = int(10000)  # Probe size threshold, in nt (Ps > 1)
    Ph = 0.1  # Maximum hole size in probe as fraction of probe size
    Po = 0.5  # Probe oligo intersection threshold for path reduction

    def _assert(self):
        ass.ert_type(self.N, int, "N")
        ass.ert_nonNeg(self.N, "N")

        ass.ert_type(self.D, int, "D")
        ass.ert_nonNeg(self.D, "D")

        ass.ert_type(self.Tr, float, "Tr")
        ass.ert_nonNeg(self.Tr, "Tr")

        ass.ert_type(self.Ps, int, "Ps")
        if self.Ps <= 1:
            raise AssertionError

        ass.ert_type(self.Ph, float, "Ph")
        ass.ert_inInterv(self.Ph, 0, 1, "Ph")

        ass.ert_type(self.Po, float, "Po")
        ass.ert_inInterv(self.Po, 0, 1, "Po")

    @staticmethod
    def get_paths(A):
        path_set = set()
        path_start_set = set()

        idxs = [
            i
            for i in range(A.shape[0] - 1)
            if i in path_start_set or len(A[i, i + 1].nonzero()[0]) != 0
        ]

        for i in idxs:
            path = [i]
            k = i
            if 1 not in A[k, (k + 1) :]:
                continue

            j = A[k, (k + 1) :].argmax() + k + 1
            path.append(j)

            while (j + 1) < A.shape[0]:
                k = j
                if 1 not in A[k, (k + 1) :]:
                    break

                j = A[k, (k + 1) :].argmax() + k + 1
                path.append(j)

            path_set.add(tuple(path))
            path_start_set.add(path[0])

        return path_set

    @staticmethod
    def get_non_overlapping_paths(oData, D):
        # Gets all paths of consecutive non-overlapping oligos with minimum
        # distance equal to D
        ass.ert_type(oData, pd.DataFrame, "oData")

        start_positions = oData["start"].values
        end_positions = oData["end"].values

        edges = [
            np.logical_or(
                end_positions + D < oData["start"].values[i],
                start_positions - D >= oData["end"].values[i],
            )
            for i in range(oData.shape[0])
        ]
        return OligoPathBuilder.get_paths(np.vstack(edges).astype("i"))

    def __size_paths(self, path_set, exit_polls):
        sized_paths = set()
        for path in set(path_set):
            if self.N > len(path):
                exit_polls["N"] = exit_polls["N"] + 1
                continue
            if self.N == len(path):
                sized_paths.add(path)
            else:
                for j in range(len(path) - self.N + 1):
                    subpath = path[j : (j + self.N)]
                    sized_paths.add(subpath)
        return sized_paths

    def filter_paths(self, path_set, oData):
        # Selects oligo paths based on length, melting temperature, size, and
        # presence of gaps.
        ass.ert_type(oData, pd.DataFrame, "oData")

        exit_polls = {"P": 0, "N": 0, "S": 0, "H": 0, "T": 0}

        selected_paths = set()
        for path in list(self.__size_paths(path_set, exit_polls)):
            if path in selected_paths:
                continue
            passed, comment = self.__path_passes(path, oData)
            exit_polls[comment] = exit_polls[comment] + 1
            if not passed:
                continue
            selected_paths.add(path)

        comment = "".join(f"{r}{c}" for (c, r) in exit_polls.items())
        return (list(selected_paths), comment)

    def __path_passes(self, path, oData):
        if not isinstance(path, list):
            path = list(path)
        pData = oData.iloc[path, :]
        probe_size = pData["end"].values[-1]
        probe_size -= pData["start"].values[0]
        if self.Ps < probe_size:
            return (False, "S")
        dData = pData["start"].values[1:] - pData["end"].values[:-1]
        max_hole_size = dData.max()
        if self.Ph * probe_size < max_hole_size:
            return (False, "H")
        Tm_range = pData["Tm"].max() - pData["Tm"].min()
        if 2 * self.Tr < Tm_range:
            return (False, "T")
        return (True, "P")

    @staticmethod
    def convert_paths_to_probes(path_list, oData):
        # Converts a list of paths into a list of probes
        ass.ert_type(path_list, list, "path list")
        probe_list = []
        if len(path_list) != 0:
            for path in path_list:
                probe_list.append(OligoProbeBuilder.path2probe(path, oData))
        return probe_list

    @staticmethod
    def path2probe(path, oData):
        # Convert an oligo path into an OligoProbe
        return OligoProbe(oData.iloc[list(path), :])


class OligoProbeBuilder(OligoPathBuilder):
    """Class to build OligoProbe objects."""

    def __init__(self):
        super(OligoProbeBuilder, self).__init__()
        self.k = None
        self.F = [0, 99]  # Threshold on number of off-targets (range)
        self.Gs = [0.0, 0.5]  # dG of SS either as kcal/mol (negative)
        #  or as fraction dG of hybrid (0<=Gs<=1)
        #  (range)
        self.Ot = 0.1  # Step for oligo score relaxation

    @property
    def config(self):
        config = cp.ConfigParser()
        config["AIM"] = {"Oligo(s) number": self.N}
        if not isinstance(self.k, type(None)):
            config["AIM"]["Oligo length (nt)"] = str(self.k)
        config["OLIGO FILTERS"] = {
            "Off-target threshold": self.F,
            "Secondary structure dG threshold": self.Gs,
            "Oligo score relaxation step": self.Ot,
        }
        config["PROBE FILTERS"] = {
            "Melting temperature range (degC)": self.Tr,
            "Min. consecutive oligo distance (nt)": self.D,
            "Probe size threshold": self.Ps,
            "Maximum hole size": self.Ph,
        }
        return config

    def _assert(self):
        OligoPathBuilder._assert(self)

        if not isinstance(self.k, type(None)):
            ass.ert_type(self.k, int, "k")
            ass.ert_nonNeg(self.k, "k")
            if (self.k + self.D) * self.N > self.Ps:
                raise AssertionError

        if isinstance(self.F, tuple):
            self.F = list(self.F)
        ass.ert_type(self.F, list, "F")
        if len(self.F) != 2:
            raise AssertionError
        for i in range(2):
            ass.ert_type(self.F[i], int, f"F[{i}]")
            if self.F[i] < 0:
                raise AssertionError
        if self.F[1] < self.F[0]:
            raise AssertionError

        if isinstance(self.Gs, tuple):
            self.Gs = list(self.Gs)
        ass.ert_type(self.Gs, list, "Gs")
        if len(self.Gs) != 2:
            raise AssertionError
        for i in range(2):
            ass.ert_type(self.Gs[i], float, f"Gs[{i}]")
            if self.Gs[i] > 1:
                raise AssertionError
        if not all(np.array(self.Gs) < 0) and not all(np.array(self.Gs) >= 0):
            raise AssertionError
        if self.Gs[0] >= 0:
            if self.Gs[1] < self.Gs[0]:
                raise AssertionError
        elif self.Gs[1] > self.Gs[0]:
            raise AssertionError

        ass.ert_type(self.Ot, float, "Ot")
        if self.Ot <= 0 or self.Ot > 1:
            raise AssertionError

    def get_prologue(self):
        s = "* OligoProbeBuilder *\n\n"
        s += f"Aim to build probes with {self.N} oligos each.\n"
        s += f"Off-target threshold range set at {self.F}.\n"
        s += "Threshold on the delta free energy of the most stable"
        if isinstance(self.Gs[0], int):
            s += f" secondary structure set at range {self.Gs} kcal/mol.\n"
        else:
            s += " secondary structure\nset at range"
            s += f" {[t*100 for t in self.Gs]}% of the delta free energy of"
            s += " hybridization.\n"

        s += f"\nMelting temperature range of {2*self.Tr} degC.\n"
        s += "Minimum distance between consecutive oligos in a probe"
        s += f" set at {self.D} nt.\n"
        s += f"Probe size threshold set at {self.Ps} nt.\n"
        s += "Reducing probes when oligo intersection fraction is equal to"
        s += f" or greater than {self.Po}.\n"

        return s

    def start(self, oGroup, window, cfr_step, logger):
        self._assert()
        return self.__build_probe_candidates(oGroup, window, cfr_step, logger)

    def __focus_oligos(self, window, oGroup):
        if np.isnan(window["cfr_start"]):
            oGroup.focus_all()
        else:
            oGroup.set_focus_window(window["cfr_start"], window["cfr_end"])
            oGroup.expand_focus_to_n_oligos(self.N)
        return oGroup

    def __build_probe_candidates(self, oGroup, window, cfr_step, logger):
        # Applies oligo filters to the oligo group,
        # expands the focus group if needed, and build probe candidates
        oGroup = self.__focus_oligos(window, oGroup)
        probe_list = self.__explore_filter(oGroup, logger)
        if not np.isnan(window["cfr_start"]):
            while len(probe_list) == 0:
                oGroup.reset_threshold()
                if oGroup.focus_window_size >= self.Ps:
                    oGroup.discard_focused_oligos_safeN(self.N - 1, self.D)
                if not oGroup.expand_focus_by_step(cfr_step):
                    break
                probe_list = self.__explore_filter(oGroup, logger)

        return probe_list

    def __build_probe_list(self, oGroup, noligos, score_thr, max_score, logger):
        nOligos_in_focus_window = oGroup.get_n_focused_oligos(True)

        probe_list = self.__get_non_overlapping_probes(
            oGroup.get_focused_oligos(True), logger
        )
        while len(probe_list) == 0:
            score_thr += self.Ot
            if score_thr > max_score:
                break

            oGroup.apply_threshold(score_thr)
            nOligosUsable = oGroup.get_n_focused_oligos(True)
            if nOligosUsable in [noligos, 0]:
                continue

            logger.info(
                " ".join(
                    [
                        f"Relaxed oligo score threshold to {score_thr:.3f}",
                        f"({nOligosUsable} oligos usable)...",
                    ]
                )
            )
            probe_list = self.__get_non_overlapping_probes(
                oGroup.get_focused_oligos(True), logger
            )

            if nOligosUsable == nOligos_in_focus_window:
                logger.warning("All oligos included. Score relaxation ineffective.")
                break
            noligos = nOligosUsable

        return probe_list

    def __explore_filter(self, oGroup, logger):
        # Explores the 0-to-max_score score threshold range and stops as soon as
        # one probe candidate passes all user-defined thresholds

        if oGroup.focus_window_size < self.Ps:
            max_score = 0
            logger.warning(
                " ".join(
                    [
                        "Score relaxation deactivated when focus region",
                        "size is smaller than probe size threshold.",
                    ]
                )
            )
        else:
            max_score = 1

        score_thr = 0
        oGroup.apply_threshold(score_thr)
        nOligos_prev_score_thr = oGroup.get_n_focused_oligos(True)
        logger.info(
            " ".join(
                [
                    f"Set oligo score threshold at {score_thr:.3f}",
                    f"({nOligos_prev_score_thr} oligos usable)...",
                ]
            )
        )

        if max_score == 0 and nOligos_prev_score_thr == 0:
            return []

        while nOligos_prev_score_thr == 0 and score_thr <= max_score - self.Ot:
            score_thr += self.Ot
            oGroup.apply_threshold(score_thr)
            nOligos_prev_score_thr = oGroup.get_n_focused_oligos(True)
            logger.info(
                " ".join(
                    [
                        f"Relaxed oligo score threshold to {score_thr:.3f}",
                        f"({nOligos_prev_score_thr} oligos usable)...",
                    ]
                )
            )

        return self.__build_probe_list(
            oGroup, nOligos_prev_score_thr, score_thr, max_score, logger
        )

    def __get_non_overlapping_probes(self, oData, logger, verbosity=True):
        # Builds non overlapping oligo paths from an oligo data table, filters
        # them based on class attributes, and then converts the remaining ones
        # to OligoProbe objects.

        paths = self.get_non_overlapping_paths(oData, self.D)
        if len(paths) == 0:
            logger.warning("No oligo paths found.")
            return []

        pathMaxLen = np.max([len(p) for p in list(paths)])
        if verbosity:
            nPaths = len(paths)
            logger.info(
                " ".join(
                    [
                        f"Found {nPaths} sets with up to {pathMaxLen}",
                        "non-overlapping oligos.",
                    ]
                )
            )
        if pathMaxLen < self.N:
            logger.warning("Longest path is shorter than requested. Skipped.")
            return []

        paths, comment = self.filter_paths(paths, oData)
        if verbosity:
            logger.info(
                " ".join(
                    [
                        f"{len(paths)}/{nPaths} oligo paths remaining",
                        f"after filtering. ({comment})",
                    ]
                )
            )

        return self.convert_paths_to_probes(paths, oData)

    def __init_probe_selection(self, sorted_probes):
        selected_probes = []
        probe_ref = sorted_probes[0]
        for probe in sorted_probes[1:-1]:
            n_shared_oligos = probe_ref.count_shared_oligos(probe)

            if probe_ref.n_oligos == n_shared_oligos:
                raise Exception("Encountered probe duplicates!")

            if self.Po * self.N <= n_shared_oligos:
                probe_ref = self.select_probe_from_pair(probe_ref, probe)
            else:
                selected_probes.append(probe_ref)
                probe_ref = probe

        if probe_ref not in selected_probes:
            selected_probes.append(probe_ref)

        return selected_probes, probe_ref

    def reduce_probe_list(self, probe_list):
        try:
            if len(probe_list) == 0:
                return []
            sorted_probes = sorted(probe_list, key=lambda p: p.range[0])

            selected_probes, probe_ref = self.__init_probe_selection(sorted_probes)

            n_shared_oligos = probe_ref.count_shared_oligos(sorted_probes[-1])
            if self.Po * self.N > n_shared_oligos:
                selected_probes.append(sorted_probes[-1])

            return selected_probes
        except Exception as e:
            print(e)
            raise

    @staticmethod
    def select_probe_from_pair(probeA, probeB):
        if probeA.size < probeB.size:
            return probeA
        if np.diff(probeA.tm_range)[0] < np.diff(probeB.tm_range)[0]:
            return probeA
        if probeA.spread / probeA.d_mean < probeB.spread / probeB.d_mean:
            return probeA
        return probeB

    @staticmethod
    def import_probes(ipath):
        # Imports output from given directory path
        if not os.path.isdir(ipath):
            raise AssertionError

        pPath = os.path.join(ipath, "probe_paths.tsv")
        if not os.path.isfile(pPath):
            return []

        oPath = os.path.join(ipath, "oligos.tsv")
        if not os.path.isfile(oPath):
            return []

        oligos = pd.read_csv(oPath, "\t", index_col=0)
        paths = pd.read_csv(pPath, "\t", index_col=0)

        return [
            OligoProbe(oligos.loc[[int(o) for o in p.split(",")]])
            for p in paths.iloc[:, 0]
        ]

    @staticmethod
    def export_probes(probe_list, opath):
        # Exports a list of probes to a directory
        if not os.path.isdir(opath):
            raise AssertionError

        probe_df = pd.concat([p.featDF for p in probe_list], ignore_index=True)
        probe_df.sort_values("start", inplace=True)
        probe_df.to_csv(os.path.join(opath, "probe_feat.tsv"), "\t")

        pd.concat([p.data for p in probe_list]).drop_duplicates().to_csv(
            os.path.join(opath, "oligos.tsv"), "\t"
        )

        probe_paths = [
            [",".join(str(x) for x in probe_list[pi].data.index.tolist())]
            for pi in range(len(probe_list))
        ]

        probe_paths = pd.DataFrame(probe_paths)
        probe_paths.columns = ["cs_oligos"]
        probe_paths.to_csv(os.path.join(opath, "probe_paths.tsv"), "\t")
