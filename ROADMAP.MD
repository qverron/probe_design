# Roadmap to Making `probe_design` user-friendly and easily repreoducible

## Orginizing Files into Folders

* Moved python and R scripts into `./src`
* Moved bash files into `./shell`

## Getting rid of executability problems

Running files like `./script.py` requires file to be executable which can require additional effort. Instead using `python ./script.py` or `bash ./command.sh` makes it simpler.

Also due to [Orginizing Files into Folders](#orginizing-files-into-folders), the names were changed as such;
    1. `python ./src/script.py`
    2. `bash ./shell/command.sh`
    3. `Rscript ./src/rcode.R`

* Changed the style in the following files;
    1.`./shell/pipeline.sh`
    2.`./shell/pipeline_excluede.sh`

## Adding pyproject.toml

A toml file is required both for creating installable modules and keeps useful information. Most importantly, it is used to keep track of required libraries.

## Adding a Licence

A licence declares how a particular code could be used. Required for making modules. Added a CC NC license that prohibits commercial use.

## Path Correction

Due to [Orginizing Files into Folders](#orginizing-files-into-folders) some **path problems** may occur. (will be a priority)

Fixed the ones in;
    1. `cycling_query.py`

## Moving code from main scope into functions

Some long chunks of the code were under `if __name__ == "__main__:"` block. While these files can be used via `python script.py` the code were not accessible outside the files. Moving those chunks to functions with the same name as their files, keeps file runnable and adds additional modularity as the functions could be accessed from other scripts.

## Renaming some python files

Some files had '-' in their names which makes it hard to use them as modules.
Renamed them into;
1. `summarize_probes_cumul.py`
2. `summarize_probes_final.py`
3. `summarize_probes.py`

## Exposing functions

Added `__init__.py` into both `./src` and into main directory which exposes functions as module.
This way;
    1. Functions inside `./src`  could be accessed from main directory
    2. Functions can be used as a part of `probe_design` library. Will be useful if we make it `pip` installable.

## Resolving circular imports in ifpd2

original had some circular import problem
for modularity and cleaner way they were resolved

## Fixing import problems in ifpd2_new

`logging.py` -> `loggingg.py`
fix conflict

## adding PATH to library

To access shell or data etc., folders inside pip installed probe_design added PATH constants into `__init__.py` . They are now accessible as `prb.PATHSHELL` considering `import probe_design as prb`.

## Making libraries pip installable

`probe_design` is now pip installable
made `ifpd2q` into standalone repo
ifpd2_new -> `ifpd2q`
removed it from `probe_design`
`ifpd2q` is now pip installable

## Updating readme

installation path changed drastically. Update it into much simpler version.

## TODO
1. ~~Move args into main scope~~
2. (maybe) make sys.argv into argparse -- Decided not to pursue
3. package other installlation requirements
4. What else??

