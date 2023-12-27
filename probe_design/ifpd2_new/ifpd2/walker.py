"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import configparser as cp
import copy
import logging
import multiprocessing as mp
import numpy as np  # type: ignore
import os
import pandas as pd  # type: ignore
from pathlib import Path
import shutil
from tqdm import tqdm  # type: ignore
from typing import Dict, List, Optional, Union

from . import asserts as ass
from .loggingg import add_log_file_handler
from .database import DataBase
from .oligo import OligoGroup
from .probe import OligoProbeBuilder
from .walker2 import ChromosomeWalker


class GenomicWindowSet:
    """Genomic window manager."""

    _window_sets = None

    C = "chr18"  # Chromosome
    S = int(3000000)  # Region start coordinate (included)
    E = int(3500000)  # Region end coordinate (excluded)

    X = 20  # Number of probes to design
    Ws: Optional[int] = None  # Window size (used when X is not provided)
    Wh: Optional[float] = 0.1  # Window shift (as a percentage of the window size)

    Rs: Union[int, float] = int(8000)  # Region focus size, either in nt (Rs > 1)
    #  or fraction of Ws (0<Rs<=1)
    #  When provided in nt, it is applied only if Rs<Ws
    Rt: Union[int, float] = int(1000)  # Region focus step, either in nt (Rt > 1)
    #  or fraction of Rs (0<Rt<=1),
    #  for focus region expansion
    _growing = False

    @property
    def window_sets(self):
        return self._window_sets

    @window_sets.setter
    def window_sets(self, window_sets):
        self._window_sets = pd.DataFrame(np.vstack(window_sets))
        self._window_sets.columns = [
            "start",
            "mid",
            "end",
            "cfr_start",
            "cfr_end",
            "w",
            "s",
        ]
        self._window_sets.sort_values("end", inplace=True)
        self._window_sets.reset_index(inplace=True, drop=True)

    @property
    def wid(self):
        return self._w

    @property
    def current_window(self):
        return self.window_sets.iloc[self.wid, :]

    @property
    def reached_last_window(self):
        return self._reached_last_window

    @property
    def growing(self):
        return self._growing

    def _assert(self):
        ass.ert_type(self.C, str, "C")
        ass.ert_type(self.S, int, "S")
        ass.ert_nonNeg(self.S, "S", True)
        ass.ert_type(self.E, int, "E")
        ass.ert_nonNeg(self.E, "E", True)
        if self.S > self.E:
            raise AssertionError

        ass.ert_multiTypes(self.X, [int, type(None)], "X")
        ass.ert_multiTypes(self.Ws, [int, type(None)], "Ws")
        if isinstance(self.X, int):
            ass.ert_type(self.Ws, type(None), "Ws")
            ass.ert_nonNeg(self.X, "X")
        else:
            ass.ert_type(self.Ws, int, "Ws")

        ass.ert_multiTypes(self.Ws, [int, type(None)], "Ws")
        if isinstance(self.Ws, int):
            ass.ert_type(self.X, type(None), "X")
            ass.ert_nonNeg(self.Ws, "Ws")
        else:
            ass.ert_type(self.X, int, "X")

        ass.ert_type(self.Wh, float, "Wh")
        ass.ert_inInterv(self.Wh, 0, 1, "Wh")

        ass.ert_multiTypes(self.Rs, [int, float], "Rs")
        if self.Rs > 1:
            self.Rs = int(self.Rs)
        else:
            ass.ert_inInterv(self.Rs, 0, 1, "Rs")

        ass.ert_multiTypes(self.Rt, [int, float], "Rt")
        if self.Rt > 1:
            self.Rt = int(self.Rt)
        else:
            ass.ert_inInterv(self.Rt, 0, 1, "Rt")

    def _init_windows(self):
        # Build windows and central focus regions (CFR)
        if not os.path.isdir(os.path.join(self.out_path, "window_sets")):
            os.mkdir(os.path.join(self.out_path, "window_sets"))

        if self.S == self.E and isinstance(self.Ws, type(None)):
            raise AssertionError(
                " ".join(
                    [
                        "During full-chromosome search, provide a window size.",
                        "I.e., it is not possible to design X probes.",
                    ]
                )
            )

        if isinstance(self.X, int):
            self.Ws = (
                self.E - self.S
                if self.X == 1
                else np.floor((self.E - self.S) / (self.X + 1)).astype("i")
            )

        self._w = 0  # Window ID
        self._reached_last_window = False

        if self.S == self.E:
            self._growing = True
            self.window_sets = self.__mk_first_window()
        else:
            self.window_sets = self.__mk_all_window_sets()

    def __update_focus(self):
        if isinstance(self.Rs, float):
            self.Rs = int(self.Rs * self.Ws)
        if isinstance(self.Rt, float):
            self.Rt = int(self.Rs * self.Rt)

    def __mk_all_window_sets(self):
        # Prepare all window sets in a region of interest
        window_starts = np.floor(np.arange(self.S, self.E, self.Ws)).astype("i")
        if (self.E - self.S) % self.Ws != 0 or len(window_starts) != 1:
            window_starts = window_starts[:-1]

        window_mids = (window_starts + self.Ws / 2).reshape((window_starts.shape[0], 1))
        window_borders = np.transpose(
            np.vstack([window_starts, window_starts + self.Ws])
        )

        nWindows = window_borders.shape[0]

        self.__update_focus()
        if self.Rs < self.Ws:
            central_regions = np.hstack(
                [
                    np.floor(window_mids - self.Rs / 2),
                    np.floor(window_mids + self.Rs / 2),
                ]
            )
        else:
            nans = np.array([np.nan for x in range(window_mids.shape[0])])
            central_regions = np.vstack([nans, nans]).transpose()

        window_sets = []
        for i in range(len(np.arange(0, 1, self.Wh))):
            winID = np.array(range(nWindows)).reshape((nWindows, 1))
            setID = np.repeat(i, nWindows).reshape((nWindows, 1))
            window_sets.append(
                np.hstack(
                    [
                        window_borders + i * self.Wh * self.Ws,
                        window_mids + i * self.Wh * self.Ws,
                        central_regions + i * self.Wh * self.Ws,
                        winID,
                        setID,
                    ]
                )[:, (0, 2, 1, 3, 4, 5, 6)]
            )

        return window_sets

    def __mk_first_window(self):
        # Initialize only the first window
        mid = self.S + self.Ws / 2
        if self.Rs < self.Ws:
            cstart, cend = (np.floor(mid - self.Rs / 2), np.floor(mid + self.Rs / 2))
        else:
            cstart, cend = (np.nan, np.nan)
        return [[self.S, mid, self.S + self.Ws, cstart, cend, 0, 0]]

    def _add_window(self):
        # Add a window, assigning it to the first set with no overlaps
        new_window = self.window_sets.iloc[-1, :].copy()
        new_window.iloc[:5] = new_window.iloc[:5] + self.Wh * self.Ws

        gData = np.array(
            [(n, g["end"].max()) for n, g in self.window_sets.groupby("s")]
        )
        non_overlapping_group_ids = np.where(gData[:, 1] <= new_window[0])[0]

        new_window_id = 0
        new_set_id = self.window_sets["s"].max() + 1
        if len(non_overlapping_group_ids) != 0:
            new_set_id = gData[non_overlapping_group_ids].min()
            new_window_id = (
                self.window_sets[self.window_sets["s"] == new_set_id]["w"].max() + 1
            )

        new_window.iloc[5] = new_window_id
        new_window.iloc[6] = new_set_id
        self.window_sets = [self.window_sets.values, new_window.values.reshape((1, 7))]

    def go_to_next_window(self):
        if self._growing:
            self._add_window()
            self._w += 1
        else:
            if self.wid < self.window_sets.shape[0] - 1:
                self._w += 1
            if self.wid >= self.window_sets.shape[0] - 1:
                self._reached_last_window = True

    def export_window_set(self):
        self.window_sets.loc[
            self.window_sets["s"] == self.current_window["s"], :
        ].to_csv(os.path.join(self.window_set_path, "windows.tsv"), "\t")


class Walker(GenomicWindowSet):
    """Walker walks through the oligos stored in an ifpd2 database,
    assigns them to windows based on user-defined parameters, and then builds
    probe candidates.
    """

    out_path = "."
    reuse = False
    _threads = 1
    __db: DataBase

    __current_oligos: List = []
    __walk_results: Dict = {}

    def __init__(self, db_path, logger=logging.getLogger()):
        GenomicWindowSet.__init__(self)
        self.__db = DataBase(db_path)
        self.logger = logger

    @property
    def threads(self):
        return self._threads

    @threads.setter
    def threads(self, threads):
        ass.ert_type(threads, int, threads)
        threads = max(1, min(threads, mp.cpu_count()))
        self._threads = threads

    @property
    def current_oligos(self):
        return self.__current_oligos

    def remove_oligos_starting_before(self, pos):
        self.__current_oligos = [o for o in self.current_oligos if o["start"] >= pos]

    @property
    def window_set_path(self):
        return os.path.join(
            self.out_path, "window_sets", f"set_{int(self.current_window['s'])}"
        )

    @property
    def window_path(self):
        return os.path.join(
            self.window_set_path, f"window_{int(self.current_window['w'])}"
        )

    @property
    def window_tag(self):
        window = self.current_window
        return f"{int(window['s'])}.{int(window['w'])}"

    @property
    def window_range(self):
        window = self.current_window
        return f"[{int(window['start'])}:{int(window['end'])}]"

    @property
    def walk_results(self):
        return self.__walk_results

    @property
    def config(self):
        config = cp.ConfigParser()

        config["AIM"] = {"Region": f"{self.C}:{self.S}-{self.E}"}
        if not isinstance(self.X, type(None)):
            config["AIM"]["Probe(s) number"] = str(self.X)

        config["WINDOWS"] = {}
        if not isinstance(self.Ws, type(None)):
            config["WINDOWS"]["Window size"] = str(self.Ws)
        config["WINDOWS"].update(
            {
                "Window step": str(self.Wh),
                "Focus region size": str(self.Rs),
                "Focus region step": str(self.Rt),
            }
        )
        return config

    def _assert(self):
        GenomicWindowSet._assert(self)
        if not os.path.isdir(self.out_path):
            raise AssertionError

    def start(self, *args, start_from_nt: int = 0, end_at_nt: int = -1, **kwargs):
        self._assert()
        self._init_windows()
        self.print_prologue()
        self.__start_walk(
            *args, start_from_nt=start_from_nt, end_at_nt=end_at_nt, **kwargs
        )

    def print_prologue(self):

        s = "* Walker *\n\n"
        s += f"Threads: {self.threads}\n"
        s += f"Database: '{self.__db.path}'\n"
        s += f"Region of interest: {self.C}:{self.S}-{self.E}\n"

        if not isinstance(self.X, type(None)):
            nWindows = int(self.window_sets["w"].max() + 1)
            s += f"Aim to scout {nWindows} windows.\n\n"

        s += f"Using a central focus region of {self.Rs} nt,"
        s += f" in windows of size {self.Ws} nt,\n"
        s += f"built with a shift of {self.Ws*self.Wh} nt ({self.Wh*100}%).\n"
        if not isinstance(self.X, type(None)):
            nSets = int(self.window_sets["s"].max() + 1)
            s += f"Thus, a total of {nSets} window sets will be explored.\n"

        self.logger.info(s)

    @property
    def walk_destination(self):
        if self.E == self.S:
            return np.inf
        return self.E

    def __finished_parsing(self, record, DBHpb, args, kwargs):
        # Return:
        #   None to stop parsing
        #   True to parse in current window
        #   Parsed output to append it

        if record["start"] < self.current_window["end"]:
            return (None, "continue")
        # DBHpb.clear()

        parsed_output = self.__process_window_async(
            self.current_oligos,
            self.current_window,
            self.fprocess,
            self.fpost,
            *args,
            opath=self.window_path,
            loggerName=self.logger.name,
            **kwargs,
        )

        if self.reached_last_window:
            self.logger.info("Reached last window")
            return (None, "break")

        self.go_to_next_window()
        self._preprocess_window()
        self._load_windows_until_next_to_do()
        if self.reached_last_window and self._window_done():
            self.logger.info("Reached last window and done")
            return (None, "break")

        if len(self.current_oligos) != 0:
            self.remove_oligos_starting_before(self.current_window["start"])

        return (parsed_output, "append")

    def __parse_line(self, record, DBHpb, args, kwargs):
        if record["start"] >= self.current_window["start"]:
            parsing_output, parsing_status = self.__finished_parsing(
                record, DBHpb, args, kwargs
            )
            if parsing_status != "continue":
                return (parsing_output, parsing_status)
            if record["end"] > self.walk_destination:  # End reached
                self.logger.info("Reached destination")
                return (None, "break")

            self.fparse(record, *args, **kwargs)

            if not np.isnan(record["score"]):
                self.current_oligos.append(record)

            self.rw += 1
        return (None, "continue")

    def __parse_database(self, DBHpb, args, kwargs):
        exec_results = []
        for record in DBHpb:
            parsing_output, parsing_status = self.__parse_line(
                record, DBHpb, args, kwargs
            )

            if parsing_status == "break":
                break
            if parsing_status == "append":
                exec_results.append(parsing_output)
            self.r += 1

        if len(self.current_oligos) != 0:
            exec_results.append(
                self.__process_window_async(
                    self.current_oligos,
                    self.current_window,
                    self.fprocess,
                    self.fpost,
                    *args,
                    opath=self.window_path,
                    loggerName=self.logger.name,
                    **kwargs,
                )
            )
        return exec_results

    def __end_walk(self, exec_results):
        for promise in exec_results:
            s, w, results = promise.get()
            if s in self.__walk_results.keys():
                self.__walk_results[s][w] = results
            else:
                self.__walk_results[s] = {w: results}

    def __start_walk(
        self, *args, start_from_nt: int = 0, end_at_nt: int = -1, **kwargs
    ):
        self.pool = mp.Pool(max(1, min(self.threads, mp.cpu_count())))
        self.logger.info(f"Prepared a pool of {self.threads} threads.")

        self._preprocess_window()
        self._load_windows_until_next_to_do()

        if self.reached_last_window and self._window_done():
            self.logger.info("All windows pre-processed. Skipped database walk.")
            return

        self.r = 0  # Record ID
        self.rw = 0  # Walk step counter

        DBHpb = tqdm(
            ChromosomeWalker(self.__db, self.C.encode()).buffer(
                start_from_nt, end_at_nt
            ),
            total=self.__db.chromosome_recordnos[self.C.encode()],
            desc="Parsing records",
            leave=False,
        )
        parsing_output = self.__parse_database(DBHpb, args, kwargs)
        self.logger.info(f"Parsed {self.rw}/{self.r} records.")

        self.__end_walk(parsing_output)
        DBHpb.close()
        self.pool.close()

    def _window_done(self):
        s = int(self.current_window["s"])
        w = int(self.current_window["w"])
        if s in self.walk_results.keys() and w in self.walk_results[s].keys():
            return isinstance(self.walk_results[s][w], list)
        return False

    def _load_windows_until_next_to_do(self):
        while self._window_done() and not self.reached_last_window:
            self.go_to_next_window()
            self._preprocess_window()

    def _preprocess_window_do(self):
        win_file_path = os.path.join(self.window_path, "window.tsv")
        win_done = os.path.isfile(os.path.join(self.window_path, ".done"))

        if os.path.isfile(win_file_path) and win_done:
            win = pd.read_csv(win_file_path, sep="\t", header=None, index_col=0)

            if (win.transpose().values[0][1:] == self.current_window.values).all():
                self.logger.info(
                    " ".join(
                        [
                            "Re-using previous results for window",
                            f"{self.window_tag} {self.window_range}",
                        ]
                    )
                )

                sid = int(self.current_window["s"])
                if sid not in self.walk_results.keys():
                    self.__walk_results[sid] = {}
                wid = int(self.current_window["w"])
                self.__walk_results[sid][wid] = self.fimport(self.window_path)
                return
        shutil.rmtree(self.window_path)
        os.mkdir(self.window_path)

    def _preprocess_window(self):
        # Preprocess current window
        # Triggers import if previously run and matching current window

        if not os.path.isdir(self.window_set_path):
            os.mkdir(self.window_set_path)
        self.export_window_set()

        if not os.path.isdir(self.window_path):
            os.mkdir(self.window_path)
        else:
            if not self.reuse:
                shutil.rmtree(self.window_path)
                os.mkdir(self.window_path)

            self._preprocess_window_do()

        self.current_window.to_csv(
            os.path.join(self.window_path, "window.tsv"), sep="\t", index=True
        )
        with open(os.path.join(self.window_path, "walker.config"), "w+") as CPH:
            self.config.write(CPH)

    @staticmethod
    def process_window_async(
        oligos,
        window,
        fprocess,
        fpost,
        *args,
        N=1,
        opath=None,
        loggerName=None,
        **kwargs,
    ):
        # Wrapper of process_window function, for parallelization
        window_tag = f"{int(window['s'])}.{int(window['w'])}"

        logger_name = f"ifpd2-window-{window_tag}"
        logger = logging.getLogger(logger_name)
        logger.propagate = False
        logger.setLevel(logging.DEBUG)
        logger.handlers = []
        add_log_file_handler("{0}/{1}.log".format(opath, "window"), logger.name)

        mainLogger = logging.getLogger(loggerName)
        mainLogger.info(f"Window {window_tag} sent to pool.")
        results = Walker.process_window(
            oligos,
            window,
            fprocess,
            fpost,
            *args,
            N=N,
            opath=opath,
            loggerName=logger.name,
            **kwargs,
        )
        mainLogger.info(f"Processor pool returned window {window_tag}.")
        return results

    def __process_window_async(self, *args, **kwargs):
        if self.pool is not None:
            return self.pool.apply_async(Walker.process_window_async, args, kwargs)
        return None

    @staticmethod
    def process_window(
        oligos,
        window,
        fprocess,
        fpost,
        *args,
        N=1,
        opath=None,
        loggerName=None,
        **kwargs,
    ):
        # Process oligos from window using fprocess. Then, post-process them
        # with fopost. Requires at least N oligos to proceed. If opath is
        # specified, a ".done" file is touched upon successful postprocessing.
        logger = logging.getLogger(loggerName)
        kwargs["loggerName"] = loggerName

        if len(oligos) >= N:
            oGroup = OligoGroup(oligos, logger)
            logger.info(
                " ".join(
                    [
                        f"Retrieved {oGroup.data.shape[0]} oligos for",
                        f"window {int(window['s'])}.{int(window['w'])}",
                        f"[{int(window['start'])}:{int(window['end'])}]",
                    ]
                )
            )
            results = fprocess(oGroup, window, *args, **kwargs)
        else:
            logger.warning(
                " ".join(
                    [
                        f"Window {int(window['s'])}.{int(window['w'])}",
                        "does not have enough oligos",
                        f"{len(oligos)}/{N}, skipped.",
                    ]
                )
            )
            results = []

        status, results = fpost(results, opath, *args, **kwargs)
        if status and opath is not None:
            Path(os.path.join(opath, ".done")).touch()

        return (int(window["s"]), int(window["w"]), results)

    @staticmethod
    def fparse(record, opb=None, *args, **kwargs):
        record.add_score(opb.F, opb.Gs)
        # if record["off_target_no"]<100:
        #     print((record._data, opb.F, opb.Gs))
        #     print(record.to_dataframe())

    @staticmethod
    def fprocess(oGroup, window, *args, **kwargs):
        opb = copy.copy(kwargs["opb"])
        if not isinstance(opb, OligoProbeBuilder):
            raise AssertionError
        logger = logging.getLogger(kwargs["loggerName"])
        probe_list = opb.start(oGroup, window, kwargs["cfr_step"], logger)
        reduced_probe_list = opb.reduce_probe_list(probe_list)
        logger.info(
            f"Reduced from {len(probe_list)} to " + f"{len(reduced_probe_list)} probes."
        )
        return reduced_probe_list

    @staticmethod
    def fimport(path, *args, **kwargs):
        return OligoProbeBuilder.import_probes(path)

    @staticmethod
    def fpost(results, opath, *args, **kwargs):
        if not isinstance(kwargs["opb"], OligoProbeBuilder):
            raise AssertionError
        logger = logging.getLogger(kwargs["loggerName"])
        if len(results) == 0:
            logger.critical(f"Built {len(results)} oligo probe candidates")
            logger.handlers = logger.handlers[:-1]
        else:
            logger.info(f"Built {len(results)} oligo probe candidates")
            logger.handlers = logger.handlers[:-1]
            OligoProbeBuilder.export_probes(results, opath)
        with open(os.path.join(opath, "builder.config"), "w+") as CPH:
            kwargs["opb"].config.write(CPH)
        return (True, results)
