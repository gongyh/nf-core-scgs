import argparse
import tempfile
from contextlib import contextmanager
from MGPGtools.utils.customHelpFormatter import CustomHelpFormatter


@contextmanager
def subparser(parser, name, desc):
    yield parser.add_parser(name, conflict_handler="resolve", help=desc, formatter_class=CustomHelpFormatter)


@contextmanager
def arg_group(parser, name):
    yield parser.add_argument_group(name)


def __db(group, required):
    group.add_argument(
        "-db",
        type=str,
        default=None,
        required=required,
        help="Path of pangenome database",
    )


def __rank(group, required):
    group.add_argument(
        "-rank",
        type=str,
        default=None,
        required=required,
        help="Rank of pangenome, choose one of 'domain', 'phylum', 'class', 'order', 'family', 'genus'",
    )


def __treeGene(group):
    group.add_argument(
        "-gene",
        type=str,
        default=None,
        help="Gene information",
    )


def __treeGenesFile(group):
    group.add_argument(
        "-genesFile",
        type=str,
        default=None,
        help="A file contains genes information",
    )


def __name(group, required):
    group.add_argument("-name", type=str, default=None, required=required, help="Name of pangenome")


def __outName(group):
    group.add_argument("-outName", type=str, default=None, help="Genus of pangenome")


def __workdir(group, required):
    group.add_argument(
        "-workdir",
        type=str,
        default=None,
        help="Directory to store tmp files",
    )


def __width(group):
    group.add_argument(
        "-x",
        "--width",
        type=int,
        default=1500,
        help="Set the width in pixels of the output image",
    )


def __height(group):
    group.add_argument(
        "-y",
        "--height",
        type=int,
        default=500,
        help="Set the height in pixels of the output image",
    )


def __path_range(group):
    group.add_argument(
        "-r",
        "--path-range",
        type=str,
        default=None,
        help="Nucleotide range to visualize: STRING=[PATH:]start-end. *-end for [0,end]; start-* for [start,pangenome_length]. If no PATH is specified, the nucleotide positions refer to the pangenome's sequence (i.e., the sequence obtained arranging all the graph's node from left to right).",
    )


def __gene(group):
    group.add_argument(
        "-gene",
        type=str,
        default=None,
        help="CDS information",
    )


def __cds(group):
    group.add_argument(
        "-cds",
        type=str,
        default=None,
        help="CDS information",
    )


def __outdir(group, required):
    group.add_argument(
        "-outdir",
        type=str,
        default=None,
        required=required,
        help="Directory to output files",
    )


def __fasta(group):
    group.add_argument(
        "-fasta",
        type=str,
        default=None,
        help="fasta file",
    )


def __coreGenes(group):
    group.add_argument(
        "-coreGenes",
        type=str,
        default=None,
        help="core genes file",
    )


def __label(group):
    group.add_argument(
        "-label",
        type=str,
        default=None,
        help="label of mapped sample",
    )


def __threads(group):
    group.add_argument("-t", "--threads", default=1, type=int, help="Number of CPUs to use")


def __debug(group):
    group.add_argument(
        "-debug",
        action="store_true",
        default=False,
        help="Create intermediate files for debugging purposes",
    )


def __help(group):
    group.add_argument("-h", "--help", action="help", help="Show help message")


def get_main_parser():
    main_parser = argparse.ArgumentParser(prog="pantools", add_help=False, conflict_handler="resolve")
    sub_parsers = main_parser.add_subparsers(help="--", dest="subparser_name")

    with subparser(sub_parsers, "stat", "static basic information of pangenome.") as parser:
        with arg_group(parser, "required named arguments") as grp:
            __db(grp, required=True)
            __outdir(grp, required=True)
            __rank(grp, required=True)
            __name(grp, required=True)
        with arg_group(parser, "optional arguments") as grp:
            __help(grp)

    with subparser(sub_parsers, "viz", "pangenome visualization") as parser:
        with arg_group(parser, "required named arguments") as grp:
            __db(grp, required=True)
            __outdir(grp, required=True)
            __name(grp, required=True)
        with arg_group(parser, "optional arguments") as grp:
            __outName(grp)
            __threads(grp)
            __help(grp)
            __width(grp)
            __height(grp)
            __path_range(grp)

    with subparser(
        sub_parsers,
        "search",
        "search pangenome information include: ref_genome, chromosome, size...",
    ) as parser:
        with arg_group(parser, "required named arguments") as grp:
            __db(grp, required=True)
            __name(grp, required=True)
        with arg_group(parser, "optional arguments") as grp:
            __gene(grp)
            __cds(grp)
            __help(grp)

    with subparser(
        sub_parsers,
        "describe",
        "describe information of gfa files",
    ) as parser:
        with arg_group(parser, "required named arguments") as grp:
            __db(grp, required=True)
            __name(grp, required=True)
            __outdir(grp, required=True)
        with arg_group(parser, "optional arguments") as grp:
            __threads(grp)
            __help(grp)

    with subparser(
        sub_parsers,
        "tree",
        "Constructing a phylogenetic tree of certain genes",
    ) as parser:
        with arg_group(parser, "required named arguments") as grp:
            __db(grp, required=True)
            __name(grp, required=True)
            __outdir(grp, required=True)
        with arg_group(parser, "optional arguments") as grp:
            __threads(grp)
            __treeGene(grp)
            __treeGenesFile(grp)
            __fasta(grp)
            __label(grp)
            __help(grp)

    with subparser(
        sub_parsers,
        "core",
        "Find Core genes of graph pangenome",
    ) as parser:
        with arg_group(parser, "required named arguments") as grp:
            __db(grp, required=True)
            __name(grp, required=True)
            __outdir(grp, required=True)
        with arg_group(parser, "optional arguments") as grp:
            __fasta(grp)
            __coreGenes(grp)
            __threads(grp)
            __help(grp)
    return main_parser
