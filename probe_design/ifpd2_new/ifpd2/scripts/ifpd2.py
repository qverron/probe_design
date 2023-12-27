"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import click  # type: ignore
from .. import __version__
from ..const import CONTEXT_SETTINGS
from ..scripts import db, extract_kmers, query, query2
import sys
import webbrowser


@click.group(
    name="ifpd2",
    context_settings=CONTEXT_SETTINGS,
    help=f"""\b
Version:    {__version__}
Author:     Gabriele Girelli
Docs:       http://ggirelli.github.io/ifpd2
Code:       http://github.com/ggirelli/ifpd2

Another iFISH probe design pipeline (II).""",
)
@click.version_option(__version__)
def main():
    """Just a hook for the entry point. Silence is golden!"""


@click.command(
    "_docs",
    help="Open online documentation on your favorite browser.",
)
def open_documentation() -> None:
    webbrowser.open("https://ggirelli.github.io/radiantkit/")
    sys.exit()


main.add_command(open_documentation,name="docs")
main.add_command(extract_kmers.main, name="extract-kmers")
main.add_command(db.run.main, name="db")
main.add_command(query.main, name="query")
main.add_command(query2.main,  name="query2")
