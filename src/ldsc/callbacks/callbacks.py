import argparse


class ChecknBlocks(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs) -> None:
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(ChecknBlocks, self).__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: int,
        option_string: str = None,
    ) -> None:
        if value <= 1:
            raise ValueError("--n-blocks must be an integer > 1.")
