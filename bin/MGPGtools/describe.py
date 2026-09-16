import logging
import os
import pandas as pd
from pathlib import Path
from MGPGtools.utils.meta import *
from MGPGtools.utils.common import run, check_directory, delete_files


class Describe(object):
    def __init__(self, options):
        self.logger = logging.getLogger("timestamp")
        self.database = options.db
        self.name = options.name
        self.outdir = options.outdir
        self.threads = 1 if options.threads is not None else options.threads
        self.meta = os.path.join(self.database, "metadata", "Metadata.tsv")
        self.gfa = os.path.join(
            self.database,
            "databases",
            "gfa",
            get_info(self.meta, self.name)["class"],
            self.name + ".gfa",
        )

    def describe_gfa(self):
        check_directory(self.outdir)
        cmd = ["gfatools", "stat", self.gfa]
        is_success, stdout, stderr = run(cmd)
        with open(os.path.join(self.outdir, self.name + ".tsv"), "w") as f:
            path_str = (
                "## databases/gfa/"
                + get_info(os.path.join(self.database, "metadata", "Metadata.tsv"), self.name)["class"]
                + "/"
                + self.name
                + ".gfa\n"
            )
            f.write(path_str)
            f.write(stdout)
