import logging
import os
import csv
from MGPGtools.exceptions import db_path_exists
from MGPGtools.utils.common import check_directory


class Stat(object):
    def __init__(self, database, rank, name, outdir):
        self.logger = logging.getLogger("timestamp")
        self.database = database
        self.rank = rank
        self.name = name
        self.outdir = outdir

    def write(self):
        check_directory(self.outdir)
        with open(
            os.path.join(self.database, "metadata", "Metadata.tsv"), "r"
        ) as f, open(os.path.join(self.outdir, self.name + ".tsv"), "w") as fout:
            reader = csv.reader(f, delimiter="\t")
            next(reader)
            header = "Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tNodes\tLinks\tNo.genomes\tRef\n"
            fout.write(header)
            for i in reader:
                if i[self.rank] == self.name:
                    fout.write("\t".join(i) + "\n")
