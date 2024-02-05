# Expose this functions to outside as a module
from .src import (
    cycling_query,
    escafish_score,
    exclude_region,
    generate_exclude,
    get_oligos,
    HUSH_feedback,
    reform_hush_combined,
    reform_hush_consec,
    reform_hush,
    select_probe,
    split_fasta,
    summarize_probes_cumul,
    summarize_probes_final,
    summarize_probes,
    download_chr_list,
    download_chr,
    download_ref_genome)


__all__ = ["cycling_query",
            "escafish_score",
            "exclude_region",
            "generate_exclude",
            "get_oligos",
            "HUSH_feedback",
            "reform_hush_combined",
            "reform_hush_consec",
            "reform_hush",
            "select_probe",
            "split_fasta",
            "summarize_probes_cumul",
            "summarize_probes_final",
            "summarize_probes",
            "download_chr_list",
            "download_chr",
            "download_ref_genome",
            ]

# CONSTANTS
import os

PATHMAIN = os.path.dirname(os.path.abspath(__file__))
PATHSRC = os.path.join(PATHMAIN, "src")
PATHSHELL = os.path.join(PATHMAIN, "shell")
PATHDATA = os.path.join(PATHMAIN, "data")
PATHNOTEBOOK = os.path.join(PATHMAIN, "notebooks")

__constants__ = ["PATHMAIN",
                "PATHSRC",
                "PATHSHELL",
                "PATHDATA",
                "PATHNOTEBOOK"]

PATHS = [PATHMAIN, PATHSRC, PATHSHELL, PATHDATA, PATHNOTEBOOK]