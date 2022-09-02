#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Command line interface to genbankpy
"""

import argparse
from pathlib import Path

from genbankpy.parser import GenBankFastaWriter


def main():
  
    parser = argparse.ArgumentParser(
        prog ='genbankpy',
        description ='Tools to download, parse and write filtered records in FASTA format out of GenBank files',
        epilog = 'Developed by Semidán Robaina Estévez (srobaina@ull.edu.es)'
    )

    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    parser._action_groups.append(optional)

    required.add_argument('--genome_ids', dest='data', type=str, required=True, nargs='+',
                        help='list of comma-separated genome ids to download')
    optional.add_argument('--outdir', dest='outdir', type=Path,
                        help='path to output directory')

    args = parser.parse_args()
    

    writer = GenBankFastaWriter.fromAccessionIDs(entry_ids=[], data_dir=args.outdir)
    writer = GenBankFastaWriter.fromSpecies(species_list=[], data_dir=args.outdir)

