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
shell_scripts_names = [x.split("/")[-1] for x in shell_scripts if "__" not in x.split("/")[-1]]
src_script_names = [x.split("/")[-1] for x in py_scripts if "__" not in x.split("/")[-1]]
notebook_scripts_names = [x.split("/")[-1] for x in notebook_scripts if "__" not in x.split("/")[-1]]

# SCRIPTS: no extension
shell_scripts_no_ext = [x.split(".")[0] for x in shell_scripts_names]
py_scripts_no_ext = [x.split(".")[0] for x in src_script_names if '.py' in x]
r_scripts_no_ext = [x.split(".")[0] for x in src_script_names if '.r' in x or '.R' in x]
notebook_scripts_no_ext = [x.split(".")[0] for x in notebook_scripts_names]

# SCRIPTS: all
ALL_SCRIPTS = shell_scripts_names + src_script_names + notebook_scripts_names
ALL_COMMANDS = shell_scripts_no_ext + py_scripts_no_ext + notebook_scripts_no_ext

def show_available() -> None:
    # ANSI escape codes for text color
    blue_code = "\033[34m"  # Dark Blue
    red_code = "\033[91m"   # Red
    green_code = "\033[92m" # Green
    reset_code = "\033[0m"  # Reset color

    print("Available scripts:\n")

    for c, s in zip(ALL_COMMANDS, ALL_SCRIPTS):
        colored_c = f"{green_code}{c}{reset_code}"
        colored_s = f"{blue_code}{s}{reset_code}"
        print(f"{colored_c:<{35}}---> ({colored_s})")

    return None


def run_command(command:str,script_arguments)->None:
    if command in shell_scripts_no_ext: # bash scripts
        subprocess.run(["bash",os.path.join(PATHSHELL, command+".sh"),*script_arguments])
    elif command in py_scripts_no_ext: # python scripts
        subprocess.run([sys.executable, os.path.join(PATHSRC, command+".py"),*script_arguments])
    elif command in notebook_scripts_no_ext: # jupyter notebooks
        subprocess.run(["jupyter", "execute", os.path.join(PATHNOTEBOOK, command+".ipynb"),*script_arguments])
    elif command in r_scripts_no_ext: # R scripts
        subprocess.run(["Rscript", os.path.join(PATHSRC, command+".r"),*script_arguments])
    else:
        print("Command not found. Please check below for the available commands.")
        show_available()
    
    return None


def main():
    if len(sys.argv) == 1:
        print("Please specify a command. Check below for the available commands.")
        print("Usage: prbdes <command> [arguments...]")
        show_available()
        return
    elif len(sys.argv) == 2 and sys.argv[1] in ["--help", "-h"]:
        print("Usage: prbdes <command> [arguments...]")
        show_available()
    else:
        command = sys.argv[1]
        script_arguments = sys.argv[2:]
        run_command(command, script_arguments)
    
    return


if __name__ == "__main__":
    main()