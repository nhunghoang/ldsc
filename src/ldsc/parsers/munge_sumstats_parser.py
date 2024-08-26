import argparse
from pathlib import Path
from .parser_interface import ParserConfig
from munge_sumstats import munge_sumstats


class MungeSumstatsconfig(ParserConfig):

    @staticmethod
    def configure_parser(
        parser: argparse.ArgumentParser, parent_parser: argparse.ArgumentParser
    ) -> argparse.ArgumentParser:
        """Add the appropriate flags and subcommands to the ldsc parser"""

        parser.add_argument(
            "--sumstats", default=None, type=Path, help="Input filename.", required=True
        )

        parser.add_argument(
            "--N",
            default=None,
            type=float,
            help="Sample size If this option is not set, will try to infer the sample "
            "size from the input file. If the input file contains a sample size "
            "column, and this flag is set, the argument to this flag has priority.",
        )
        parser.add_argument(
            "--N-cas",
            default=None,
            type=float,
            help="Number of cases. If this option is not set, will try to infer the number "
            "of cases from the input file. If the input file contains a number of cases "
            "column, and this flag is set, the argument to this flag has priority.",
        )
        parser.add_argument(
            "--N-con",
            default=None,
            type=float,
            help="Number of controls. If this option is not set, will try to infer the number "
            "of controls from the input file. If the input file contains a number of controls "
            "column, and this flag is set, the argument to this flag has priority.",
        )

        parser.add_argument(
            "--info-min", default=0.9, type=float, help="Minimum INFO score."
        )
        parser.add_argument("--maf-min", default=0.01, type=float, help="Minimum MAF.")

        # daner and --daner-n are not compatible so we are going to add them to a mutually exclusive group
        daner_group = parser.add_mutually_exclusive_group(required=False)

        daner_group.add_argument(
            "--daner",
            default=False,
            action="store_true",
            help="Use this flag to parse Stephan Ripke's daner* file format.",
        )
        daner_group.add_argument(
            "--daner-n",
            default=False,
            action="store_true",
            help="Use this flag to parse more recent daner* formatted files, which "
            "include sample size column 'Nca' and 'Nco'.",
        )

        # --no-alleles and --merge-alleles are also not compatible so we are going
        # to add them to a group
        alleles_group = parser.add_mutually_exclusive_group(required=False)

        alleles_group.add_argument(
            "--no-alleles",
            default=False,
            action="store_true",
            help="Don't require alleles. Useful if only unsigned summary statistics are available "
            "and the goal is h2 / partitioned h2 estimation rather than rg estimation.",
        )
        alleles_group.add_argument(
            "--merge-alleles",
            default=None,
            type=str,
            help="Same as --merge, except the file should have three columns: SNP, A1, A2, "
            "and all alleles will be matched to the --merge-alleles file alleles.",
        )
        parser.add_argument(
            "--n-min",
            default=None,
            type=float,
            help="Minimum N (sample size). Default is (90th percentile N) / 2.",
        )
        parser.add_argument("--chunksize", default=5e6, type=int, help="Chunksize.")

        # optional args to specify column names
        parser.add_argument(
            "--snp",
            default=None,
            type=str,
            help="Name of SNP column (if not a name that ldsc understands). NB: case insensitive.",
        )
        parser.add_argument(
            "--N-col",
            default=None,
            type=str,
            help="Name of N column (if not a name that ldsc understands). NB: case insensitive.",
        )
        parser.add_argument(
            "--N-cas-col",
            default=None,
            type=str,
            help="Name of N column (if not a name that ldsc understands). NB: case insensitive.",
        )
        parser.add_argument(
            "--N-con-col",
            default=None,
            type=str,
            help="Name of N column (if not a name that ldsc understands). NB: case insensitive.",
        )
        parser.add_argument(
            "--a1",
            default=None,
            type=str,
            help="Name of A1 column (if not a name that ldsc understands). NB: case insensitive.",
        )
        parser.add_argument(
            "--a2",
            default=None,
            type=str,
            help="Name of A2 column (if not a name that ldsc understands). NB: case insensitive.",
        )
        parser.add_argument(
            "--p",
            default=None,
            type=str,
            help="Name of p-value column (if not a name that ldsc understands). NB: case insensitive.",
        )
        parser.add_argument(
            "--frq",
            default=None,
            type=str,
            help="Name of FRQ or MAF column (if not a name that ldsc understands). NB: case insensitive.",
        )
        parser.add_argument(
            "--signed-sumstats",
            default=None,
            type=str,
            help="Name of signed sumstat column, comma null value (e.g., Z,0 or OR,1). NB: case insensitive.",
        )
        parser.add_argument(
            "--info",
            default=None,
            type=str,
            help="Name of INFO column (if not a name that ldsc understands). NB: case insensitive.",
        )
        parser.add_argument(
            "--info-list",
            default=None,
            type=str,
            help="Comma-separated list of INFO columns. Will filter on the mean. NB: case insensitive.",
        )
        parser.add_argument(
            "--nstudy",
            default=None,
            type=str,
            help="Name of NSTUDY column (if not a name that ldsc understands). NB: case insensitive.",
        )
        parser.add_argument(
            "--nstudy-min",
            default=None,
            type=float,
            help="Minimum # of studies. Default is to remove everything below the max, unless there is an N column,"
            " in which case do nothing.",
        )
        parser.add_argument(
            "--ignore",
            default=None,
            type=str,
            help="Comma-separated list of column names to ignore.",
        )
        parser.add_argument(
            "--a1-inc",
            default=False,
            action="store_true",
            help="A1 is the increasing allele.",
        )
        parser.add_argument(
            "--keep-maf",
            default=False,
            action="store_true",
            help="Keep the MAF column (if one exists).",
        )

        parser.set_defaults(func=munge_sumstats)

        return parser
