# Expose the following functions to the main package
from .cycling_query_new import cycling_query
from .escafish_score_new import escafish_score
from .exclude_region import exclude_region
from .generate_exclude import generate_exclude
from .get_oligos import get_oligos
from .HUSH_feedback import HUSH_feedback
from .reform_hush_combined_new import reform_hush_combined
from .reform_hush_consec import reform_hush_consec
from .reform_hush import reform_hush
from .select_probe import select_probe
from .split_fasta import split_fasta
from .summarize_probes_cumul import summarize_probes_cumul
from .summarize_probes_final import summarize_probes_final
from .summarize_probes import summarize_probes
from .download_chromosomes import download_chr_list, download_chr, download_ref_genome


__all__ = ["cycling_query_new",
            "escafish_score_new",
            "exclude_region",
            "generate_exclude",
            "get_oligos",
            "HUSH_feedback",
            "reform_hush_combined_new",
            "reform_hush_consec",
            "reform_hush",
            "select_probe",
            "split_fasta",
            "summarize_probes_cumul",
            "summarize_probes_final",
            "summarize_probes",
            "download_chr_list",
            "download_chr",
            "download_ref_genome"]

import os
# PATHMAIN is different from main init file
PATHMAIN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
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