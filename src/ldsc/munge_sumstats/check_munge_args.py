from typing import Any


def check_munge_args(args: list[Any]) -> None:
    if args.out is None:
        raise ValueError("No output was specified")
