#!/usr/bin/env python3

import os
import logging
from Bio import SeqIO
import argparse

def windowScheme(length, size, spacing):
    currentPosition = 0
    outputIndices = []
    while currentPosition <= length:
        outputIndices.append((currentPosition, currentPosition + size if currentPosition + size < length else length))
        currentPosition += spacing + size
        if currentPosition >= length:
            break
    return outputIndices

def chunkString(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))


logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()
## Script to output a kmer breakdown shell script to run and partition data files for use

parser = argparse.ArgumentParser(description='Tool for plotting kmer information criteria across .')
parser.add_argument('-i', '--input', dest='input_sequence', type=str, help='Input sequence file', required=True)
parser.add_argument('-f', '--format', dest='sequence_format', type=str, help='Input file type', default='fasta')
parser.add_argument('-o', '--output', dest='output', type=str, help='Output prefix', default='kmer_output')
parser.add_argument('-w', '--window-size', dest='size', type=int, help='Window size for sequences', default=1000000)
parser.add_argument('-s', '--window-overlap', dest='spacing', type=int, help='Window overlap for sequences', default=0)
parser.add_argument('-k', '--kmer-range', dest='kmers', type=str, help='kmer size scheme as: min,max,steps', default="15,25,5")
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)

args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.INFO)

if not os.path.exists(args.input_sequence):
    logger.error(f"Couldn't find input sequence file '{args.input_sequence}. Exiting.'")
    exit()

sequences = SeqIO.to_dict(SeqIO.parse(args.input_sequence, args.sequence_format))
for sequence in sequences.keys():
    schemes = windowScheme(len(sequences[sequence].seq), args.size, args.spacing)
    for scheme in schemes:
        with open(f"{args.output}_window_{scheme[0]}-{scheme[1]}.fna", "w") as fileOut:
            fileOut.write(f">{sequence}_window_{scheme[0]}-{scheme[1]}\n")
            chunks = chunkString(str(sequences[sequence].seq[scheme[0]:scheme[1]]), length=80)
            for chunk in chunks:
                fileOut.write(f"{chunk}\n")
