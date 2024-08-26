from datetime import datetime
import traceback
from ldsc.parsers import generate_parser
from ldsc.logger import LDSCLogger


def main() -> None:
    """main function that will be called for ldsc"""

    main_parser = generate_parser()

    args = main_parser.parse_args()

    start_time = datetime.now()

    # create a logger and configure it
    logger = LDSCLogger.create_logger()

    log_stem, log_prefix = args.log_filename.split(".")

    args.log_filename = f"{log_stem}_{args.func.__name__}.{log_prefix}"

    logger.configure(
        args.out.parent, args.log_filename, args.verbose, args.log_to_console
    )

    logger.print_header(args, main_parser)

    logger.info(f"Analysis started at {start_time}")
    try:
        # Now we are going to call the appropriate function for the selected subcommand
        args.func(args)

        end_time = datetime.now()

        logger.info(f"Analysis finished at {end_time}")
    except Exception as e:
        logger.critical("Encountered an error during the analysis")
        logger.critical(traceback.format_exc())
    # print(main_parser)


if __name__ == "__main__":
    main()
