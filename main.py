from datetime import datetime
import parsers  


def main() -> None:
    """main function that will be called for ldsc"""

    main_parser = parsers.generate_parser()

    args = main_parser.parse_args()

    start_time = datetime.now()

    # Now we are going to call the appropriate function for the selected subcommand
    args.func(args)

    end_time = datetime.now()

    # print(main_parser)

    

if __name__ == "__main__":
    main()
