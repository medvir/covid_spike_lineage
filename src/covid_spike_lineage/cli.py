#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 07:38:42 2021

"""
import argparse
import sys
from pkg_resources import (get_distribution, DistributionNotFound)
from argparse import RawTextHelpFormatter
import os

try:
    __version__ = get_distribution('covid_spike_lineage').version
except:
    __version__ = 'unknown'

parser = argparse.ArgumentParser(description='Detect interesting lineages from'
                                 ' a part of Spike Sanger sequences smaller than'
                                 ' (aa 403 to aa 760).'
                                     , formatter_class=RawTextHelpFormatter)

sanger_path = os.path.abspath(os.getcwd())
parser.add_argument('-d','--directory', dest='seq_dir_path', type=str, 
                    default=sanger_path,
                    help='Enter the Absolute sequence directory path')

parser.add_argument('-e','--extention', dest='ext', type=str, 
                    default='ab1',
                    #default='',
                    #required=True,
                    help='Enter the file type to be analysed (ab1 or fasta). ')

parser.add_argument('-v','--version', action='version', version=__version__)

#args = parser.parse_args()

def main(args=None):

    import logging
    import logging.handlers

    
    args = parser.parse_args(args=args)
    
    log_format = '%(levelname)s %(asctime)s %(filename)s: %(funcName)s() %(lineno)d: \t%(message)s'

    log_full_path = os.path.join(args.seq_dir_path,'S_gene_lineage.log')
    logging.basicConfig(filename=log_full_path, level=logging.INFO,
                        format=log_format, datefmt='%Y/%m/%d %H:%M:%S')
    logging.info(' '.join(sys.argv))

    from covid_spike_lineage import lineage_detection
    lineage_detection.main(args.seq_dir_path, args.ext)