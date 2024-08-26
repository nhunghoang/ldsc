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

from ldsc.logger import LDSCLogger

logger: logging.Logger = LDSCLogger.get_logger(__name__)


def check_args(args: list[Any]) -> None:
    # pull the arguments from the parser into a dictionary so then i can safely check to see which arguments are present
    function_args = vars(args)

    if function_args.get("bfile") is not None and function_args.get("l2") is None:
        raise argparse.ArgumentError("Must specify --l2 with --bfile.")

    if (
        function_args.get("cts_bin") is not None
        and function_args.get("extract") is not None
    ):
        raise argparse.ArgumentError(
            "--cts-bin and --extract are currently incompatible."
        )
    if (
        function_args.get("annot") is not None
        and function_args.get("cts_bin") is not None
    ):
        raise argparse.ArgumentError(
            "--annot and --cts-bin are currently incompatible."
        )
    if (function_args.get("cts_bin") is not None) != (
        function_args.get("cts_breaks") is not None
    ):
        raise argparse.ArgumentError(
            "Must set both or neither of --cts-bin and --cts-breaks."
        )

    if function_args.get("per_allele"):
        args.pq_exp = 1

    if function_args.get("h2") is not None and function_args.get("rg") is not None:
        raise argparse.ArgumentError("Cannot set both --h2 and --rg.")

    if (function_args.get("samp_prev") is not None) != (
        function_args.get("pop_prev") is not None
    ):
        raise argparse.ArgumentError(
            "Must set both or neither of --samp-prev and --pop-prev."
        )

    if not function_args.get("overlap_annot") or function_args.get("not_M_5_50"):
        if (
            function_args.get("frqfile") is not None
            or function_args.get("frqfile_chr") is not None
        ):
            logger.info("The frequency file is unnecessary and is being ignored.")
            args.frqfile = None
            args.frqfile_chr = None
    if function_args.get("overlap_annot") and not function_args.get("not_M_5_50"):
        if not (
            (function_args.get("frqfile") and function_args.get("ref_ld"))
            or (function_args.get("frqfile_chr") and function_args.get("ref_ld_chr"))
        ):
            raise argparse.ArgumentError(
                "Must set either --frqfile and --ref-ld or --frqfile-chr and --ref-ld-chr"
            )
