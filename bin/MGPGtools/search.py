import os
import logging
from prettytable import PrettyTable
from MGPGtools.utils.meta import *


class Search(object):
    def __init__(self, options):
        self.logger = logging.getLogger("timestamp")
        self.database = options.db
        self.name = options.name
        if options.gene is not None or options.cds is not None:
            self.region_type = "gene" if options.gene is not None else "CDS"
            self.region_name = options.gene if options.gene is not None else options.cds
        else:
            self.region_type = ""
            self.region_name = ""

    def ref_info(self):
        info, chrom_size_info = get_chromosome(self.database, self.name)
        path = info["ref"].replace(".", "#")
        tb_row_count = 1
        if len(self.region_name) != 0:
            gff_path = os.path.join(
                self.database,
                "databases",
                "gff",
                info["class"],
                info["ref"] + ".gff",
            )
            character = self.region_type.lower() + "-" + self.region_name
            id_character = "ID=" + character
            start = ""
            end = ""
            tchrom = ""
            with open(gff_path, "r") as f:
                lines = f.read().strip().split("\n")
                for line in lines:
                    if line[0] == "#":
                        continue
                    row = line.strip().split("\t")
                    if (
                        row[2] == self.region_type
                        and row[8].split(";")[0] == id_character
                    ):
                        tchrom = row[0]
                        start = row[3]
                        end = row[4]
            if len(start) == 0:
                print(
                    "{} {} not found in {}".format(
                        self.region_type, self.region_name, info["ref"]
                    )
                )
            else:
                tb = PrettyTable(
                    [
                        "name",
                        "ref",
                        self.region_type,
                        "chromosome",
                        self.region_type + "_path",
                        "start",
                        "end",
                    ]
                )
                tb.add_row(
                    [
                        self.name,
                        info["ref"],
                        self.region_name,
                        tchrom,
                        path + "#" + tchrom,
                        start,
                        end,
                    ]
                )
                print(tb)
        else:
            tb = PrettyTable(["name", "ref", "chromosomes", "size", "path"])
            for chrom, size in chrom_size_info.items():
                if tb_row_count == 1:
                    tb.add_row(
                        [self.name, info["ref"], chrom, size, path + "#" + chrom]
                    )
                else:
                    tb.add_row(["", "", chrom, size, path + "#" + chrom])
                tb_row_count += 1
            print(tb)
