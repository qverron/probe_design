import sys
import os
import argparse
import subprocess

# PATH constants
PATHMAIN = os.path.dirname(__file__)
PATHSRC = os.path.join(PATHMAIN, "src")
PATHSHELL = os.path.join(PATHMAIN, "shell")
PATHDATA = os.path.join(PATHMAIN, "data")
PATHNOTEBOOK = os.path.join(PATHMAIN, "notebooks")

# SCRIPTS: absolute paths
shell_scripts = [os.path.join(PATHSHELL, x) for x in os.listdir(PATHSHELL)]
py_scripts = [os.path.join(PATHSRC, x) for x in os.listdir(PATHSRC)]
notebook_scripts = [os.path.join(PATHNOTEBOOK, x) for x in os.listdir(PATHNOTEBOOK)]

# SCRIPTS: names
shell_scripts_names = [x.split("/")[-1] for x in shell_scripts]
py_scripts_names = [x.split("/")[-1] for x in py_scripts]
notebook_scripts_names = [x.split("/")[-1] for x in notebook_scripts]

# SCRIPTS: no extension
shell_scripts_no_ext = [x.split(".")[0] for x in shell_scripts_names]
py_scripts_no_ext = [x.split(".")[0] for x in py_scripts_names]
notebook_scripts_no_ext = [x.split(".")[0] for x in notebook_scripts_names]




def main():
    print("Hello World: SUCCESS!")
    print(shell_scripts_names)


if __name__ == "__main__":
    main()