import argparse
import sys

# add argparser that takes -l and -m and multiply them
argparser = argparse.ArgumentParser()
argparser.add_argument("-l", type=int, help="length")
argparser.add_argument("-x", type=int, help="height")
argparser.add_argument("-w", type=int, help="width")
argparser.add_argument("-s", type=str, help="shape")
args = argparser.parse_args()


if __name__ == "__main__":

    volume = args.l * args.x * args.w
    shape = args.s
    print(f"Volume of {shape} is {volume}")

    # prb test -l 10 -h 20 -w 30 -s box