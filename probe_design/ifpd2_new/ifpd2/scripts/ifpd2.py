"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import click  # type: ignore
from ifpd2 import __version__
from ifpd2.const import CONTEXT_SETTINGS
from ifpd2.scripts import db, extract_kmers, query, query2
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


main.add_command(open_documentation)
main.add_command(extract_kmers.main)
main.add_command(db.run.main)
main.add_command(query.main)
main.add_command(query2.main)
