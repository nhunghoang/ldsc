from datetime import datetime
from parsers  import generate_parser
from logger import LDSCLogger


def main() -> None:
    """main function that will be called for ldsc"""

    main_parser = generate_parser()

    args = main_parser.parse_args()

    start_time = datetime.now()

    # create a logger and configure it
    logger = LDSCLogger.create_logger()

    logger.configure(
        args.out.parent, args.log_filename, args.verbose, args.log_to_console
    )
    
    logger.print_header(args, main_parser)

    logger.info(f"Analysis started at {start_time}")

    # Now we are going to call the appropriate function for the selected subcommand
    args.func(args)

    end_time = datetime.now()

    logger.info(f"Analysis finished at {end_time}")

    # print(main_parser)

    

if __name__ == "__main__":
    main()
