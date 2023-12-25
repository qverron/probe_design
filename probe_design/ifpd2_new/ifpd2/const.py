"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

dtype_melting = {
    "name": "|S",
    "Tm_dG": "<f4",
    "Tm_dH": "<f4",
    "Tm_dS": "<f4",
    "Tm": "<f4",
    "sequence": "|S",
}
dtype_hush = {"sequence": "|S", "off_target_no": "<u8"}
dtype_secondary = {"ss_dG": "<f4"}
dtype_sequence_features = {"sequence": "|S", "gc_content": "<f4"}
dtype_header_features = {"name": "|S", "chromosome": "|S", "start": "<u4", "end": "<u4"}

DEFAULT_DATABASE_INDEX_BIN_SIZE = 100000

database_columns = [
    "name",
    "chromosome",
    "start",
    "end",
    "sequence",
    "gc_content",
    "off_target_no",
    "Tm_dG",
    "Tm_dH",
    "Tm_dS",
    "Tm",
    "ss_dG",
]
