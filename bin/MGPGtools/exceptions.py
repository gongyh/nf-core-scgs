import os
from MGPGtools.utils.exceptions import *
from MGPGtools.utils.taxonomy import Taxonomy
import MGPGtools.utils.meta as umeta


def db_path_exists(db_path):
    if not os.path.exists(db_path):
        raise incorrectDBPath("incorrect database path {})".format(db_path))


def rankName_exists(db, rank, rankName):
    rankList = umeta.rank_list(os.path.join(db, "metadata", "Metadata.tsv"), rank)
    if rankName not in rankList:
        raise pangenomeNotFound("No {} pangenome is found".format(rankName))


def stat_args_judgment(args):
    if args.rank not in Taxonomy.rank_labels:
        raise rankNotFound(
            "{} is not the rank ('domain','phylum','class','order','family','genus')".format(
                args.rank
            )
        )
    rankName_exists(args.db, args.rank, args.name)


def viz_args_judgment(args):
    rankName_exists(args.db, "genus", args.name)
