#!/usr/bin/env python
"""
(c) 2014 Brendan Bulik-Sullivan and Hilary Finucane

LDSC is a command line tool for estimating
    1. LD Score
    2. heritability / partitioned heritability
    3. genetic covariance / correlation

"""

import logging
from typing import Any

import argparse

from logger import LDSCLogger

logger: logging.Logger = LDSCLogger.get_logger(__name__)


def check_args(args: list[Any]) -> None:
    if args.bfile is not None and args.l2 is None:
            raise argparse.ArgumentError("Must specify --l2 with --bfile.")

    if args.cts_bin is not None and args.extract is not None:
        raise argparse.ArgumentError("--cts-bin and --extract are currently incompatible.")
    if args.annot is not None and args.cts_bin is not None:
        raise argparse.ArgumentError("--annot and --cts-bin are currently incompatible.")
    if (args.cts_bin is not None) != (args.cts_breaks is not None):
        raise argparse.ArgumentError(
            "Must set both or neither of --cts-bin and --cts-breaks."
        )
    
    if args.per_allele:
        args.pq_exp = 1

    if args.h2 is not None and args.rg is not None:
        raise argparse.ArgumentError("Cannot set both --h2 and --rg.")
    

    if (args.samp_prev is not None) != (args.pop_prev is not None):
        raise argparse.ArgumentError(
            "Must set both or neither of --samp-prev and --pop-prev."
        )

    if not args.overlap_annot or args.not_M_5_50:
        if args.frqfile is not None or args.frqfile_chr is not None:
            logger.info("The frequency file is unnecessary and is being ignored.")
            args.frqfile = None
            args.frqfile_chr = None
    if args.overlap_annot and not args.not_M_5_50:
        if not (
            (args.frqfile and args.ref_ld)
            or (args.frqfile_chr and args.ref_ld_chr)
        ):
            raise argparse.ArgumentError(
                "Must set either --frqfile and --ref-ld or --frqfile-chr and --ref-ld-chr"
            )
