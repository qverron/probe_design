import argparse
import sys

# takes a command with undefined number of arguments
def print_args(*many)->None:

    for one in many:
        print(one)
    
    return None

if __name__ == "__main__":

    print(sys.argv)
    print_args(*sys.argv)

