#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Command line interface to genbankpy
"""

import os
import sys
import argparse
from filtersam.filtersam import filterSAM


def main():
  
    parser = argparse.ArgumentParser(
        prog ='genbankpy',
        description ='Tools to download, parse and write FASTA out of GenBank files',
        epilog = 'Developed by Semidán Robaina Estévez (srobaina@ull.edu.es)'
    )

    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    parser._action_groups.append(optional)

    required.add_argument('--hmm', dest='hmm', type=str, required=True,
                        help='path to tigrfam or pfam model')
    required.add_argument('--in', dest='data', type=str, required=True,
                        help='path to peptide database')
    optional.add_argument('--outdir', dest='outdir', type=str,
                        help='path to output directory')
    optional.add_argument('--prefix', dest='prefix', type=str,
                        default='',
                        help='prefix to be added to output files')

    args = parser.parse_args()
    bam = os.path.abspath(args.bam)

    if not os.path.isfile(bam):
        print('Specified bam file does not exist')
        sys.exit()

    if args.identity is not None:
        print(f'Filtering by percent identity at {args.identity}%')
        filterSAM(input_path=bam, output_path=args.out, filter_by='identity',
                  cutoff=args.identity, n_processes=args.processes)
    if args.matched is not None:
        print(f'Filtering by percent of matched sequence at {args.matched}%')
        filterSAM(input_path=bam, output_path=args.out, filter_by='matched',
                  cutoff=args.matched, n_processes=args.processes)
    if args.identity is None and args.matched is None:
        print(f'Defaulting to filtering by percent identity at 95%')
        filterSAM(input_path=bam, output_path=args.out, filter_by='identity',
                  cutoff=95.0, n_processes=args.processes)  