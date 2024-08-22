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

from log import CustomLogger

logger: logging.Logger = CustomLogger.get_logger(__name__)


__version__ = "2.0.2"
MASTHEAD = "*********************************************************************\n"
MASTHEAD += "* LD Score Regression (LDSC)\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane\n"
MASTHEAD += "* Broad Institute of MIT and Harvard / MIT Department of Mathematics\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "*********************************************************************\n"



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





#     """main function that will have the logic for ldsc"""
# args = parser.parse_args()

# log = Logger(args.out.name +  ".log")

# log2 = CustomLogger.create_logger()

# log2.configure(args.out.parent, f"{args.out.name}.log2", 1, True)

# defaults = vars(parser.parse_args(""))
# opts = vars(args)
# non_defaults = [x for x in list(opts.keys()) if opts[x] != defaults[x]]
# header = MASTHEAD
# header += "Call: \n"
# header += "./ldsc.py \\\n"
# options = [
# "--" + x.replace("_", "-") + " " + str(opts[x]) + " \\"
# for x in non_defaults
# ]
# header += "\n".join(options).replace("True", "").replace("False", "")
# header = header[0:-1] + "\n"
# log.log(header)
# log.log("Beginning analysis at {T}".format(T=time.ctime()))
# start_time = time.time()




# if __name__ == "__main__":

#     args = parser.parse_args()

#     log = Logger(args.out.name +  ".log")

#     log2 = CustomLogger.create_logger()

#     log2.configure(args.out.parent, f"{args.out.name}.log2", 1, True)

#     try:
#         defaults = vars(parser.parse_args(""))
#         opts = vars(args)
#         non_defaults = [x for x in list(opts.keys()) if opts[x] != defaults[x]]
#         header = MASTHEAD
#         header += "Call: \n"
#         header += "./ldsc.py \\\n"
#         options = [
#             "--" + x.replace("_", "-") + " " + str(opts[x]) + " \\"
#             for x in non_defaults
#         ]
#         header += "\n".join(options).replace("True", "").replace("False", "")
#         header = header[0:-1] + "\n"
#         log.log(header)
#         log.log("Beginning analysis at {T}".format(T=time.ctime()))
#         start_time = time.time()

#             # bad flags
#         else:
#             print(header)
#             print("Error: no analysis selected.")
#             print("ldsc.py -h describes options.")
#     except Exception:
#         ex_type, ex, tb = sys.exc_info()
#         log.log(traceback.format_exc(ex))
#         raise
#     finally:
#         log.log("Analysis finished at {T}".format(T=time.ctime()))
#         time_elapsed = round(time.time() - start_time, 2)
#         # log.log("Total time elapsed: {T}".format(T=sec_to_str(time_elapsed)))
