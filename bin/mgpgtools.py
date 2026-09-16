#!/usr/bin/env python3
import logging
import sys
from MGPGtools.cli import get_main_parser
from MGPGtools.analysis import Analysis
from MGPGtools.info import __version__, __author__


def print_help():
    print(
        """\

        ...::: MGPGtools v%s:::...

    Methods:
    describe -> basic information of gfa files
    stat -> statistic information of pangenome
    viz -> pangenome visualization
    search -> search pangenome information include: ref_genome, chromosome, size...

    Use: MGPGtools <command> -h for command specific help
    """
        % __version__
    )


def main():
    args = None
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {"-v", "--v", "-version", "--version"}:
        print(f"MGPGtools: version {__version__} {__author__}")
    elif sys.argv[1] in {"-h", "--h", "-help", "--help"}:
        print_help()
        sys.exit(0)
    else:
        args = get_main_parser().parse_args()

    analysis = Analysis()
    if args.subparser_name == "stat":
        analysis.statistic(args)
    elif args.subparser_name == "viz":
        analysis.viz(args)
    elif args.subparser_name == "search":
        analysis.search(args)
    elif args.subparser_name == "describe":
        analysis.describe(args)
    elif args.subparser_name == "tree":
        analysis.tree(args)
    elif args.subparser_name == "core":
        analysis.core(args)


if __name__ == "__main__":
    main()
