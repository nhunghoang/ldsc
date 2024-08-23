from typing import Any
from log import CustomLogger
from .default_header import MASTHEAD

class LDSCLogger(CustomLogger):

    def __init__(self) -> None:
        super().__init__()
    
    def print_header(self, args: list[Any], parser) -> None:
        self.info(MASTHEAD)

        defaults = vars(parser.parse_args(""))
        opts = vars(args)
        non_defaults = [x for x in list(opts.keys()) if opts[x] != defaults[x]]
        header = MASTHEAD
        header += "Call: \n"
        header += "./ldsc.py \\\n"
        options = [
            "--" + x.replace("_", "-") + " " + str(opts[x]) + " \\"
            for x in non_defaults
        ]

        header += "\n".join(options).replace("True", "").replace("False", "")
        header = header[0:-1] + "\n"
        self.header(header)


        