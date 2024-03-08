#!/bin/bash

mkdir data
mkdir data/rois
mkdir data/ref
mkdir HUSH


usage() {
  echo "Usage: $0 -g genome_path -r regions_path"
  echo "  -g, --genome    Path to reference genome"
  echo "  -r, --regions   Path to regions of interest"
  echo "  -h, --help      Display this help message"
  echo "Example: $0 -g /path/to/genome.fa -r /path/to/regions.tsv"
  exit 1
}


# add --genome and -g for taking path to reference genome and --regions and -r for taking path to regions of interest

while (( "$#" )); do
  case "$1" in
    -g|--genome)
      genome="$2"
      shift 2
      ;;
    -r|--regions)
      regions="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done


# create symbloic links to the reference genome and regions of interest in the data directory i n ref and rois subdirectories

ln -s $genome data/ref/genome.fa
ln -s $regions data/rois/all_regions.tsv
