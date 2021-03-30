import argparse
import sys

from metagraph.cli import workflows


def main(args=tuple(sys.argv[1:])):
    parser = argparse.ArgumentParser(description='metagraph utils')

    subparsers = parser.add_subparsers(help="Available subcommands", required=True,
                                       dest="command")

    build_parser = subparsers.add_parser("build", help="Create index")
    workflows.setup_build_parser(build_parser)

    parsed_arguments = parser.parse_args(args)

    if parsed_arguments.func:
        parsed_arguments.func(parsed_arguments)
    else:
        sys.exit("Unknown function call")


if __name__ == "__main__":
    main()
