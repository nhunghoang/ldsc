import logging
from typing import Any
from log import CustomLogger
from .default_header import MASTHEAD


class LDSCLogger(CustomLogger):

    def __init__(self, name:str, level: int = logging.NOTSET) -> None:
        super().__init__(name, level)

    def print_header(self, args: list[Any], parser) -> None:

        opts = vars(args)

        non_defaults = [x for x in list(opts.keys()) if opts[x] is not None]
        header = MASTHEAD
        header += "Call: \n"
        header += "./ldsc.py \\\n"
        options = [
            "--" + x.replace("_", "-") + " " + str(opts[x]) + " \\"
            for x in non_defaults
        ]

        header += "\n".join(options).replace("True", "").replace("False", "")
        header = header[0:-1] + "\n"
        self.info(header)


logging.setLoggerClass(LDSCLogger)    