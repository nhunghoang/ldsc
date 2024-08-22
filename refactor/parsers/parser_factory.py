import argparse
from pathlib import Path
from .parser_interface import ParserConfig
from .ldsc_parser import LDSCconfig
from .munge_sumstats_parser import MungeSumstatsconfig
from .make_annot_parser import MakeAnnotconfig
from rich_argparse import RichHelpFormatter


def _configure_parsers(parsers: list[argparse.ArgumentParser]) -> None:
    """configure the parser to have the appropriate commandline arguments based on the choices the user selects"""

    configuration_opts: dict[str, ParserConfig] = {
        "ldsc": LDSCconfig,
        "munge_sumstats": MungeSumstatsconfig,
        "make_annot": MakeAnnotconfig
    }

    for parser in parsers:
        if (config_obj := configuration_opts.get(parser.description)):
            config_obj.configure_parser(parser)
        else:
            raise ValueError("There was an issue configuring the parsers")


def generate_parser() -> argparse.ArgumentParser:
    # Creating the main parser that will be used in the code
    __version__ = "2.2.0"
    # Here we are going to create the main overall parser and 
    # the common option parser. We will configure the subparser 
    # in different scripts
    main_parser = argparse.ArgumentParser(
        description="`ldsc` is a command line tool for estimating heritability and genetic correlation from GWAS summary statistics. `ldsc` also computes LD Scores.",
        formatter_class=RichHelpFormatter,
        add_help=True,
    )

    main_parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s: {__version__}",
    )

    # We are also going to create a parser for common options
    common_parser = argparse.ArgumentParser(
        formatter_class=RichHelpFormatter, add_help=False
    )

    common_parser.add_argument(
        "-o",
        "--out",
        default="ldsc",
        type=Path,
        help="output filename prefix. If only a filename is given then it will be written to the current directory. If a full path with the output prefix is given then the file will be written there. (default: %(default)s)"
    )

    # Creating a group for logging parameters. These arguments
    # don't change how downstream functions run
    logging_group = common_parser.add_argument_group(
        "logging",
        description="parameters that affect what runtime information is printed to a file or the console",
    )

    logging_group.add_argument(
        "--verbose",
        "-v",
        default=0,
        help="verbose flag indicating if the user wants more information",
        action="count",
    )

    logging_group.add_argument(
        "--log-to-console",
        default=False,
        help="Optional flag to log to only the console or also a file",
        action="store_true",
    )

    logging_group.add_argument(
        "--log-filename",
        default="test.log",
        type=str,
        help="Name for the log output file. (default: %(default)s)",
    )

    # now we are going to instantiate the subparsers
    subparser = main_parser.add_subparsers(
        title="Commands",
        description="three different pipelines that the user can select. These pipelines are ldsc, munge_sumstats, and make_annot'.",
    )

    #each subparser has a description that will be used to determine how to properly configure the parser
    ldsc_parser = subparser.add_parser(
        name="ldsc",
        help="pipeline to estimate LD Score, heritability/partitioned heritablility, and genetic covariance",
        formatter_class=RichHelpFormatter,
        parents=[common_parser],
        description="ldsc"
    )

    munge_sumstats_parser = subparser.add_parser(
        "munge_sumstats",
        help="pipeline used to format the summary statistics file to be input for ldsc.",
        formatter_class=RichHelpFormatter,
        parents=[common_parser],
        description="munge_sumstats"
    )

    make_annot_parser = subparser.add_parser(
        "make_annot",
        help="pipeline used to make the annotation file. This pipeline requires that bedtools be installed and be a part of the system path.",
        formatter_class=RichHelpFormatter,
        parents=[common_parser],
        description="make_annot"
    )

    
    _configure_parsers([ldsc_parser, munge_sumstats_parser, make_annot_parser])

    return main_parser
    # configure_parser