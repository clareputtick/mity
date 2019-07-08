"""mity API and command-line interface"""
import argparse
import logging
import os
import sys

from . import (call, normalise)
from ._version import __version__

__all__ = []
def public(fn):
    __all__.append(fn.__name__)
    return fn

AP = argparse.ArgumentParser(
        description="mity, a Mitochondrial DNA variant analysis toolkit.",
        epilog=f"See the online manual for details: http://github.com/KCCG/mity. Version {__version__}")
AP_subparsers = AP.add_subparsers(
        help="Sub-commands (use with -h for more info)")


# call ------------------------------------------------------------------------

do_call = public(call.do_call)

def _cmd_call(args):
  """Call mitochondrial variants"""
  logging.info("mity %s", __version__)
  logging.info("Calling mitochondrial variants")

  call.do_call(args.bam, args.reference, args.prefix, args.min_mq, args.min_bq, args.min_af, args.min_ac, args.ploidy, args.normalise)

P_call = AP_subparsers.add_parser('call', help=_cmd_call.__doc__)
P_call.add_argument('bam', action = 'append', nargs='+', help = 'BAM files to run the analysis on.')
P_call.add_argument('--reference', action='store', required=True)
P_call.add_argument('--prefix', action='store', help = 'Output files will be named with PREFIX')
P_call.add_argument('--min-mapping-quality', action='store', type = int, default = 30, help = 'Exclude alignments from analysis if they have a mapping quality less than MIN_MAPPING_QUALITY. Default: 30', dest="min_mq")
P_call.add_argument('--min-base-quality', action='store', type = int, default = 20, help = 'Exclude alleles from analysis if their supporting base quality is less than MIN_BASE_QUALITY. Default: 20', dest="min_bq")
P_call.add_argument('--min-alternate-fraction', action='store', type = float, default = 0.5, help = 'Require at least MIN_ALTERNATE_FRACTION observations supporting an alternate allele within a single individual in the in order to evaluate the position. Default: 0.0001, range = [0,1]', dest="min_af")
P_call.add_argument('--min-alternate-count', action='store', type = int, default = 4, help = 'Require at least MIN_ALTERNATE_COUNT observations supporting an alternate allele within a single individual in order to evaluate the position. Default: 4', dest="min_ac")
P_call.add_argument('--ploidy', action='store', type = int, default = 2, help = 'Expected ploidy of the sample. Default: 2.')
P_call.add_argument('--normalise', action='store_true', help = 'Normalise the resulting VCF? This is currently broken')
# parser.add_argument('--parallel', action='store_true', help = 'Run freebayes in parallel.')
# parser.add_argument('--ncpu', action='store', type = int, help = 'Number of CPUs to use when running in parallel.')

P_call.set_defaults(func=_cmd_call)


# normalise ------------------------------------------------------------------------

do_normalise = public(normalise.do_normalise)

def _cmd_normalise(args):
    """Normalise & FILTER mitochondrial variants"""
    logging.info("mity %s", __version__)
    logging.info("Normalising and FILTERing mitochondrial vcf.gz file")
    
    normalise.do_normalise(args.vcf, args.outfile)

P_normalise = AP_subparsers.add_parser('normalise', help=_cmd_normalise.__doc__)
P_normalise.add_argument('--vcf', action='store', required=True, help="vcf.gz file from running mity")
P_normalise.add_argument('--outfile', action='store', required=True, help="output VCF file in bgzip compressed format")
P_normalise.set_defaults(func=_cmd_normalise)


# version ------------------------------------------------------------------------

def print_version(_args):
    """Display this program's version."""
    print(__version__)

P_version = AP_subparsers.add_parser('version', help=print_version.__doc__)
P_version.set_defaults(func=print_version)

def parse_args(args=None):
    """Parse the command line."""
    return AP.parse_args(args=args)
