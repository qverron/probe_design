import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("-l", type=int, help="length", default=1)
argparser.add_argument("-x", type=int, help="height", default=2)
argparser.add_argument("-w", type=int, help="width", default=3)
argparser.add_argument("-s", type=str, help="shape", default="box")
args = argparser.parse_args()

def main():

    green_code = "\033[92m" # Green
    reset_code = "\033[0m"  # Reset color
    red_code = "\033[91m"   # Red

    volume = args.l * args.x * args.w
    shape = args.s
    print(f"{red_code}Length: {args.l}, Height: {args.x}, Width: {args.w}{reset_code}")
    print(f"Volume of {shape} is {volume}")

    print(f"{green_code}TEST SUCCESFUL!{reset_code}")
    # prb test -l 10 -x 20 -w 30 -s box
if __name__ == "__main__":

    main()