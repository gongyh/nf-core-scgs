import logging
import os
import subprocess
from MGPGtools.utils.meta import *
from MGPGtools.utils.exceptions import pathNotFound
from MGPGtools.exceptions import db_path_exists
from MGPGtools.utils.common import run, check_directory, delete_files


class Odgi(object):
    def __init__(self, database, name, outdir, outName=None, threads=1):
        """Initalization

        Args:
            name (str): Genus name
            outdir (str): Output directory
            outName (str, optional): Output name. Defaults to None.
            threads (int, optional): Number of threads to use. Defaults to 1.
        """
        self.logger = logging.getLogger("timestamp")
        self.database = database
        self.name = name
        self.outdir = outdir
        if outName is None:
            self.outName = name
        else:
            self.outName = outName
        self.threads = threads
        # check_on_path("odgi")

    def build(self):
        check_directory(self.outdir)
        gfa_path = os.path.join(
            self.database,
            "databases",
            "gfa",
            get_info(os.path.join(self.database, "metadata", "Metadata.tsv"), self.name)["class"],
            self.name + ".gfa",
        )
        og_path = os.path.join(self.outdir, self.outName + ".og")
        cmd = (
            ["odgi", "build", "-g", gfa_path, "-o", og_path]
            if self.threads == 1
            else [
                "odgi",
                "build",
                "-g",
                gfa_path,
                "-o",
                og_path,
                "-t",
                str(self.threads),
            ]
        )
        is_success, stdout, stderr = run(cmd)

    def sort(self):
        og_path = os.path.join(self.outdir, self.outName + ".og")
        sort_og_path = os.path.join(self.outdir, self.outName + ".sorted.og")
        cmd = (
            ["odgi", "sort", "-i", og_path, "-P", "-Y", "-O", "-o", sort_og_path]
            if self.threads == 1
            else [
                "odgi",
                "sort",
                "-i",
                og_path,
                "-P",
                "-Y",
                "-O",
                "-o",
                sort_og_path,
                "-t",
                str(self.threads),
            ]
        )
        is_success, stdout, stderr = run(cmd)

    def viz(self, options):
        self.build()
        self.sort()
        sort_og_path = os.path.join(self.outdir, self.outName + ".sorted.og")
        png_path = os.path.join(self.outdir, self.outName + ".png")
        cmd = ["odgi", "viz", "-i", sort_og_path, "-o", png_path]
        if self.threads != 1:
            cmd.extend(["-t", str(self.threads)])
        if options.width is not None:
            cmd.extend(["-x", str(options.width)])
        if options.height is not None:
            cmd.extend(["-y", str(options.height)])
        if options.path_range is not None:
            info, chrom_size = get_chromosome(self.database, self.name)
            path_list = []
            for i in chrom_size:
                path_list.append(info["ref"].replace(".", "#") + "#" + i)
            if ":" in options.path_range:
                if options.path_range.split(":")[0] not in path_list:
                    raise pathNotFound("Path not found in {}".format(self.name))
                else:
                    cmd.extend(["-r", options.path_range])
            else:
                if options.path_range not in path_list:
                    raise pathNotFound("Path not found in {}".format(self.name))
                else:
                    path_size = "{}:0-{}".format(
                        options.path_range,
                        chrom_size[options.path_range.split("#")[-1]],
                    )
                    cmd.extend(["-r", path_size])
        is_success, stdout, stderr = run(cmd)
        delete_files(os.path.join(self.outdir, self.outName + ".sorted.og"))
        delete_files(os.path.join(self.outdir, self.outName + ".og"))
