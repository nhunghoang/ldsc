import argparse
from typing import Protocol


class ParserConfig(Protocol):

    @staticmethod
    def configure_parser(parser: argparse.ArgumentParser, parent_parser: argparse.ArgumentParser) -> None:
        ...