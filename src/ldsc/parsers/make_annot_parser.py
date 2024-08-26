import argparse
from pathlib import Path
from .parser_interface import ParserConfig
from annotate import make_annot_files


class MakeAnnotconfig(ParserConfig):

    @staticmethod
    def configure_parser(
        parser: argparse.ArgumentParser, parent_parser: argparse.ArgumentParser = None
    ) -> argparse.ArgumentParser:
        """Add the appropriate flags and subcommands to the ldsc parser"""

        parser.add_argument(
            "--gene-set-file",
            type=Path,
            help="a file of gene names, one line per gene.",
        )
        parser.add_argument(
            "--gene-coord-file",
            type=Path,
            default="ENSG_coord.txt",
            help="a file with columns GENE, CHR, START, and END, where START and END are base pair coordinates of TSS and TES. This file can contain more genes than are in the gene set. We provide ENSG_coord.txt as a default.",
        )
        parser.add_argument(
            "--windowsize",
            type=int,
            help="how many base pairs to add around the transcribed region to make the annotation?",
        )
        parser.add_argument(
            "--bed-file",
            type=Path,
            help="the UCSC bed file with the regions that make up your annotation",
        )
        parser.add_argument(
            "--nomerge",
            action="store_true",
            default=False,
            help="don't merge the bed file; make an annot file wi    th values proportional to the number of intervals in the bedfile overlapping the SNP.",
        )
        parser.add_argument(
            "--bimfile",
            type=Path,
            help="plink bim file for the dataset you will use to compute LD scores.",
        )
        parser.add_argument(
            "--annot-file", type=Path, help="the name of the annot file to output."
        )

        parser.set_defaults(func=make_annot_files)

        return parser
