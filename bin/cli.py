#!/usr/bin/env python3

import time
from enum import Enum
from pathlib import Path
import typer
from Bio import SeqIO
import subprocess
import numpy as np


def real_split(
    fa,
    ann,
    level_Bacteria,
    level_Eukaryota,
    bac_out_dir,
    euk_out_dir,
    gff=None,
    ko_file=None,
):
    """
    Split by annotation
    """
    record_dict = SeqIO.to_dict(SeqIO.parse(fa, "fasta"))

    ctg_genes = {}
    if gff is None:  # skip splitting gene ids
        pass
    else:  # parse the gff file
        with open(gff) as gfh:
            for line in gfh:
                cl = line.strip().split("\t")
                if len(cl) >= 9:  # gene annotation line
                    ctg_id = cl[0]
                    anno8 = cl[8].split(";")[0]
                    if anno8.startswith("ID=") and anno8.endswith("_gene"):  # correct line
                        gid = anno8[3 : len(anno8) - 5]
                        if ctg_id in ctg_genes.keys():  # add the new id
                            ctg_genes[ctg_id].append(gid)
                        else:
                            ctg_genes[ctg_id] = [gid]

    gene_ko = {}
    if ko_file is None:
        pass
    else:
        with open(ko_file) as kfh:
            for line in kfh:
                cl = line.strip().split("\t")
                if len(cl) == 2:  # correct annotation
                    if cl[0] not in gene_ko.keys():
                        gene_ko[cl[0]] = [cl[1]]
                    else:
                        gene_ko[cl[0]].append(cl[1])

    annCol = -1  # which colum corresponds to the specified level
    eukAnnCol = -1  # which colum corresponds to the superkindom level
    eukCol = -1
    with open(ann) as fh:
        for line in fh:
            if line[0:2] == "##":  # comment line, skip
                continue
            elif line[0:6] == "# name":  # header line
                cl = line.strip().split("\t")
                for item in cl:
                    il = item.split(".")
                    if len(il) == 3 and il[0] == level_Bacteria and il[1] == "t":  # eg. order.t.12
                        annCol = int(il[2].rstrip("%s")) - 1
                    if len(il) == 3 and il[0] == level_Eukaryota and il[1] == "t":
                        eukAnnCol = int(il[2].rstrip("%s")) - 1
                    if len(il) == 3 and il[0] == "superkingdom" and il[1] == "t":
                        eukCol = int(il[2].rstrip("%s")) - 1
                continue
            cl = line.strip().split("\t")
            contig_id = cl[0]
            ctg_id_short = contig_id.split("_length_")[0]
            superkingdom = cl[eukCol]
            if superkingdom == "Eukaryota":
                eukAnnotation = cl[eukAnnCol].replace("/", "_")
                if "_" in eukAnnotation:
                    eukAnnotation = eukAnnotation.split("_")[0]
                euk_gname = eukAnnotation.replace(" ", "_") + ".gids"
                euk_kname = eukAnnotation.replace(" ", "_") + ".ko"
                euk_gname_path = euk_out_dir.joinpath(euk_gname)
                euk_kofh = euk_out_dir.joinpath(euk_kname).open("a")
                with euk_gname_path.open("a") as f2:
                    if ctg_id_short in ctg_genes.keys():
                        for g in ctg_genes[ctg_id_short]:
                            f2.write(g + "\n")
                            if g in gene_ko.keys():
                                for v in gene_ko[g]:
                                    euk_kofh.write(g + "\t" + v + "\n")
                euk_kofh.close()
            else:
                bacAnnotation = cl[annCol].replace("/", "_")
                if "_" in bacAnnotation:
                    bacAnnotation = bacAnnotation.split("_")[0]
                bac_gname = bacAnnotation.replace(" ", "_") + ".gids"
                bac_kname = bacAnnotation.replace(" ", "_") + ".ko"
                bac_kofh = bac_out_dir.joinpath(bac_kname).open("a")
                bac_gname_path = bac_out_dir.joinpath(bac_gname)
                with bac_gname_path.open("a") as f1:
                    if ctg_id_short in ctg_genes.keys():
                        for g in ctg_genes[ctg_id_short]:
                            f1.write(g + "\n")
                            if g in gene_ko.keys():
                                for v in gene_ko[g]:
                                    bac_kofh.write(g + "\t" + v + "\n")
                bac_kofh.close()
            assert annCol > 0
            assert eukAnnCol > 0
            assert eukCol > 0
            bac_fname = bacAnnotation.replace(" ", "_") + ".fasta"
            euk_fname = eukAnnotation.replace(" ", "_") + ".fasta"
            if superkingdom == "Eukaryota":
                fname_path = euk_out_dir.joinpath(euk_fname)
            else:
                fname_path = bac_out_dir.joinpath(bac_fname)
            with fname_path.open("a") as f:
                SeqIO.write(record_dict[contig_id], f, "fasta")


class TaxaLevels(str, Enum):
    superkingdom = "superkingdom"
    phylum = "phylum"
    order = "order"
    family = "family"
    genus = "genus"
    species = "species"


state = {"verbose": False}
__version__ = "1.0"
APP_NAME = "scgs-cli"

app = typer.Typer()
tools_app = typer.Typer()
app.add_typer(tools_app, name="tools", help="Misc tools.")


@tools_app.command("scgs_split", help="Split assembly genome according to taxa")
def tools_split(
    results_dir: Path = typer.Option(
        "./results/",
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=False,
        readable=True,
        resolve_path=True,
        show_default=True,
    ),
    level_Bacteria: TaxaLevels = typer.Option(TaxaLevels.genus, case_sensitive=False, show_default=True),
    level_Eukaryota: TaxaLevels = typer.Option(TaxaLevels.genus, case_sensitive=False, show_default=True),
    output_dir: Path = typer.Option(
        "./split/",
        exists=False,
        file_okay=False,
        dir_okay=True,
        writable=True,
        readable=True,
        resolve_path=True,
        show_default=True,
    ),
    force: bool = typer.Option(False, "--force", "-f"),
):
    typer.echo(f"Checking the existence of blob, prokka and spades results.")
    samples = []

    blob_dir = results_dir.joinpath("blob")
    if not blob_dir.is_dir():
        typer.secho(
            f"Taxa annotations not found, please check {blob_dir} .",
            fg=typer.colors.RED,
        )
        raise typer.Abort()
    else:  # check subdir
        for child in blob_dir.iterdir():
            if child.is_dir():  # sample
                samples.append(child.parts[-1])
    typer.echo(f"INFO: {len(samples)} samples detected.")

    spades_dir = results_dir.joinpath("spades")
    if not spades_dir.is_dir():
        typer.secho(
            f"Spades assemblies not found, please check {spades_dir} .",
            fg=typer.colors.RED,
        )
        raise typer.Abort()
    else:  # check the existence of all genome assemblies
        for sample in samples:
            ass = spades_dir.joinpath(sample + ".ctg200.fasta")
            if not ass.exists():
                typer.secho(f"genome assembly file {ass} not found.", fg=typer.colors.RED)
                raise typer.Abort()

    prokka_dir = results_dir.joinpath("prokka")
    anno_exist = False
    if prokka_dir.exists() and prokka_dir.is_dir():  # exist the dir, ok
        typer.echo(f"INFO: gene annotations found, will also split gene ids.")
        anno_exist = True
    else:
        typer.echo(f"INFO: gene annotations not found, will skip splitting gene ids.")
        anno_exist = False

    kofam_dir = results_dir.joinpath("kofam")
    ko_exist = False
    if kofam_dir.exists() and kofam_dir.is_dir():
        typer.echo(f"INFO: KO annotations found, will also split gene kos.")
        ko_exist = True
    else:
        typer.echo(f"INFO: KO annotations not found, will skip splitting gene kos.")
        ko_exist = False

    typer.secho(f"Done!", fg=typer.colors.GREEN)

    # try to create outdir
    try:
        if force and output_dir.exists() and output_dir.is_dir():
            output_dir.rename(str(output_dir) + "%d" % (int(time.time())))
            output_dir.mkdir(exist_ok=False)
        else:
            output_dir.mkdir(exist_ok=True)
    except:
        typer.secho(f"Can not create output directory.", fg=typer.colors.RED)
        raise typer.Abort()

    # process split one by one
    with typer.progressbar(samples, label="Processing") as progress:
        for sample in progress:
            # typer.echo(sample)
            fa = spades_dir.joinpath(sample + ".ctg200.fasta")
            out_bac_subdir = output_dir.joinpath(sample + "_" + level_Bacteria + "_Bacteria")
            out_bac_subdir.mkdir(exist_ok=True)  # create subdir to store Bacteria fastas
            out_euk_subdir = output_dir.joinpath(sample + "_" + level_Eukaryota + "_Eukaryota")
            out_euk_subdir.mkdir(exist_ok=True)  # create subdir to store Eukaryota fastas
            blob_sub = blob_dir.joinpath(sample)
            blob_table = None
            for child in blob_sub.iterdir():
                if child.match(sample + ".blobDB*.table.txt"):
                    blob_table = child
            if blob_table is None:
                typer.secho(
                    f"Can not find annotation table for sample {sample}.",
                    fg=typer.colors.RED,
                )
                raise typer.Abort()
            gff = None
            sample_anno = prokka_dir.joinpath(sample)
            if anno_exist and sample_anno.exists() and sample_anno.is_dir():
                sample_gff = sample_anno.joinpath(sample + ".gff")
                if sample_gff.exists() and sample_gff.is_file():  # find annotation gff file
                    gff = sample_gff
            ko_file = None
            if ko_exist:
                ko_file = kofam_dir.joinpath(sample + "_KOs_mapper.txt")
            # perform split
            real_split(
                fa,
                blob_table,
                level_Bacteria,
                level_Eukaryota,
                out_bac_subdir,
                out_euk_subdir,
                gff,
                ko_file,
            )
        typer.secho(f"\nFinished.", fg=typer.colors.GREEN)


@tools_app.command("scgs_checkm", help="Perform CheckM for (splitted) draft assemblies")
def tools_checkm(
    fastas_dir: Path = typer.Option(
        ...,
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=False,
        readable=True,
        resolve_path=True,
        show_default=True,
    ),
    suffix: str = typer.Option("fasta", show_default=True),
    genus: str = typer.Option("lineage_wf", show_default=True),
    output_dir: Path = typer.Option(
        "./checkm_out/",
        exists=False,
        file_okay=False,
        dir_okay=True,
        writable=True,
        readable=True,
        resolve_path=True,
        show_default=True,
    ),
    threads: int = typer.Option(16, show_default=True),
    checkm_table: Path = typer.Option(
        "./checkm.txt",
        exists=False,
        file_okay=True,
        dir_okay=False,
        writable=True,
        readable=True,
        resolve_path=True,
        show_default=True,
    ),
    force: bool = typer.Option(False, "--force", "-f"),
):
    typer.echo(f"Checking output directory.")

    if output_dir.exists() and output_dir.is_dir():
        if force:
            output_dir.rename(str(output_dir) + "%d" % (int(time.time())))
            typer.secho(f"Rename checkm output dir.", fg=typer.colors.RED)
        else:
            typer.secho(
                f"Output directory already exist, please move/delete and try again.",
                fg=typer.colors.RED,
            )
            raise typer.Abort()

    if checkm_table.exists() and checkm_table.is_file():
        if force:
            checkm_table.rename(str(checkm_table) + "%d" % (int(time.time())))
            typer.secho(f"Rename checkm result file.", fg=typer.colors.RED)
        else:
            typer.secho(f"Output checkm results already exist.", fg=typer.colors.RED)
            raise typer.Abort()

    subprocess.check_call(
        " ".join(
            [
                str(Path(__file__).resolve().parent.joinpath("checkm.sh")),
                str(fastas_dir),
                str(suffix),
                str(genus),
                str(output_dir),
                str(threads),
                str(checkm_table),
            ]
        ),
        shell=True,
    )
    typer.secho(f"Finished", fg=typer.colors.GREEN)


@tools_app.command("scgs_fastANI", help="Perform FastANI for genomes")
def tools_fastANI(
    query_genome: Path = typer.Option(
        None,
        "--query",
        "-q",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        show_default=False,
        help="Query genome (fasta/fastq)[.gz]",
    ),
    query_list: Path = typer.Option(
        None,
        "--queryList",
        "--ql",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        show_default=False,
        help="A file containing list of query genomes, one genome per line. * Incompatible with: --query",
    ),
    ref_genome: Path = typer.Option(
        None,
        "--ref",
        "-r",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        show_default=False,
        help="Reference genome (fasta/fastq)[.gz]",
    ),
    ref_list: Path = typer.Option(
        None,
        "--refList",
        "--rl",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        show_default=False,
        help="A file containing list of reference genomes, one genome per line. * Incompatible with: --ref",
    ),
    output_file: Path = typer.Option(
        "fastani.out",
        "--output",
        "-o",
        exists=False,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        show_default=True,
        help="output file name",
    ),
    threads: int = typer.Option(
        1,
        "--threads",
        "-t",
        show_default=True,
        help="Thread count for parallel execution.",
    ),
    visualize: bool = typer.Option(False, "--visualize", help="Output mappings and visualization."),
):
    """
    This is a wrapper script for fastANI. FastANI is a fast alignment-free implementation
    for computing whole-genome Average Nucleotide Identity (ANI) between genomes.
    """

    if query_genome is not None and query_list is not None:
        typer.secho(
            f"Error: only one of --query or --queryList can be set.",
            fg=typer.colors.RED,
        )
        raise typer.Abort()

    if ref_genome is not None and ref_list is not None:
        typer.secho(f"Error: only one of --ref or --refList can be set.", fg=typer.colors.RED)
        raise typer.Abort()

    if visualize and (query_genome is None or ref_genome is None):
        typer.secho(
            f"Error: visualize can only be enabled for one to one genome comparison",
            fg=typer.colors.RED,
        )
        raise typer.Abort()

    params_str = ""

    if query_genome is not None:
        params_str += "--query " + str(query_genome)
    elif query_list is not None:
        params_str += "--queryList " + str(query_list)
    else:
        typer.secho(
            f"Error: query genome(s) need to be set by --query or --queryList.",
            fg=typer.colors.RED,
        )
        raise typer.Abort()

    if ref_genome is not None:
        params_str += " --ref " + str(ref_genome)
    elif ref_list is not None:
        params_str += " --refList " + str(ref_list)
    else:
        typer.secho(
            f"Error: reference genome(s) need to be set by --ref or --refList.",
            fg=typer.colors.RED,
        )
        raise typer.Abort()

    if visualize:
        params_str += " --visualize"

    subprocess.check_call(
        " ".join(
            [
                "fastANI",
                params_str,
                "--output",
                str(output_file),
                "--threads",
                str(threads),
            ]
        ),
        shell=True,
    )
    if visualize:
        visfile = Path(str(output_file) + ".visual")
        if visfile.exists() and visfile.is_file():  # can be used to draw plot
            subprocess.check_call(
                " ".join(
                    [
                        str(Path(__file__).resolve().parent.joinpath("fastANI_vis.R")),
                        str(query_genome),
                        str(ref_genome),
                        str(visfile),
                    ]
                ),
                shell=True,
            )
        else:  # error
            typer.secho(
                f"Error: mapping file (with .visual extension) is not generated.",
                fg=typer.colors.RED,
            )

    typer.secho(f"\nFinished.", fg=typer.colors.GREEN)


@tools_app.command("scgs_roary", help="Perform Roary for genomes")
def tools_roary(
    input_dir: Path = typer.Option(
        ...,
        "--input",
        "-i",
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=False,
        readable=True,
        resolve_path=True,
        show_default=False,
        help="Directory with gff files.",
    ),
    threads: int = typer.Option(8, "--threads", "-t", show_default=True, help="Number of threads."),
    kraken_db: Path = typer.Option(
        None,
        "--kraken_db",
        "-k",
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=False,
        readable=True,
        resolve_path=True,
        show_default=False,
        help="Path to Kraken database for QC.",
    ),
    output_dir: Path = typer.Option(
        ".",
        "--output",
        "-f",
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=True,
        readable=True,
        resolve_path=True,
        show_default=True,
        help="Output directory.",
    ),
):
    """
    This is a wrapper script for Roary. All GFF3 files created by Prokka are valid with Roary
    and this is the recommended way of generating the input files.
    """
    typer.echo(f"Checking input.")
    gffs = []

    for child in input_dir.iterdir():
        if child.is_file() and child.suffix == ".gff":
            gffs.append(child)

    if len(gffs) < 3:
        typer.secho(f"Error: at least 3 gff files needed.", fg=typer.colors.RED)
        raise typer.Abort()

    params_str = "-r -e --mafft"
    if kraken_db is not None:
        params_str += " -qc -k " + str(kraken_db)

    typer.echo(f"Run roary.")
    subprocess.check_call(
        " ".join(
            [
                "roary",
                params_str,
                "-f",
                str(output_dir),
                "-p",
                str(threads),
                str(input_dir) + "/*.gff",
            ]
        ),
        shell=True,
    )
    gpa = output_dir.joinpath("gene_presence_absence.csv")
    if not gpa.exists():
        typer.secho(f"Error: gene_presence_absence.csv not generated", fg=typer.colors.RED)
        raise typer.Abort()
    aln = output_dir.joinpath("core_gene_alignment.aln")
    if not aln.exists():
        typer.secho(f"Error: no core gene alignment file generated.", fg=typer.colors.RED)
        raise typer.Abort()

    typer.echo(f"Generate a newick tree.")
    subprocess.check_call(
        " ".join(["fasttree", "-nt", "-gtr", str(aln), ">", str(output_dir) + "/tree.newick"]),
        shell=True,
    )

    tree = output_dir.joinpath("tree.newick")
    if not tree.exists():
        typer.secho(f"Error: tree file not generated.", fg=typer.colors.RED)
        raise typer.Abort()

    typer.echo(f"Plot figures.")
    subprocess.check_call(
        " ".join(
            [
                str(Path(__file__).resolve().parent.joinpath("roary_plots.py")),
                str(tree),
                str(gpa),
            ]
        ),
        shell=True,
    )

    typer.secho(f"Finished", fg=typer.colors.GREEN)


@tools_app.command("scgs_scoary", help="Perform Scoary for genomes")
def tools_scoary(
    traits: Path = typer.Option(
        ...,
        "--traits",
        "-t",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        show_default=False,
        help="Input trait table (csv), with trait presence: 1, absence: 0.",
    ),
    gene_presence: Path = typer.Option(
        ...,
        "--genes",
        "-g",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        show_default=False,
        help="Input gene presence/absence table (csv) from ROARY.",
    ),
    correction: str = typer.Option(
        "I",
        "--correction",
        "-c",
        show_default=True,
        help="Apply the indicated filtration measures (I B BH PW EPW or P), separate multiple measures by space, e.g. EPW P",
    ),
    pvalue_cutoff: str = typer.Option(
        "0.05",
        "--p_value_cutoff",
        "-p",
        show_default=True,
        help="P-value cut-off / alpha level, can be multiple values for different correction method, separated by space.",
    ),
    permute: int = typer.Option(
        None,
        "--permute",
        "-e",
        show_default=False,
        help="Perform N number of permutations of the significant results post-analysis. (E.g. 10000)",
    ),
    restrict_to: str = typer.Option(
        None,
        "--restrict_to",
        "-r",
        show_default=False,
        help="Use if you only want to analyze a subset of your strains. (E.g. Strain1,Strain2,Strain3)",
    ),
    threads: int = typer.Option(1, "--threads", "-t", show_default=True, help="Number of threads."),
    newicktree: Path = typer.Option(
        None,
        "--newicktree",
        "-n",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        show_default=False,
        help="Supply a custom tree (Newick format) for phylogenetic analyses instead instead of calculating it internally.",
    ),
    collapse: bool = typer.Option(False, "--collapse", help="Collapse correlated genes into merged units."),
    outdir: Path = typer.Option(
        ".",
        "--outdir",
        "-o",
        exists=False,
        file_okay=False,
        dir_okay=True,
        writable=True,
        readable=True,
        resolve_path=True,
        show_default=True,
        help="Directory to place output files.",
    ),
):
    """
    This is a wrapper script for Scoary.
    """
    typer.echo(f"Checking input.")
    cmd = "scoary -g " + str(gene_presence) + " -t " + str(traits)

    corrects = correction.strip().split(" ")
    pvalues = pvalue_cutoff.strip().split(" ")
    if len(corrects) > 0 and len(pvalues) > 0 and len(corrects) == len(pvalues):  # correct
        cmd += " -c " + corrects + " -p " + pvalues
    else:
        typer.secho(
            f"Error: check --correction and --p_value_cutoff parameters.",
            fg=typer.colors.RED,
        )
        raise typer.Abort()

    if permute is not None:
        cmd += " -e " + permute

    if restrict_to is not None:
        cmd += " -r " + restrict_to + " -w"

    cmd += " --threads " + threads

    if newicktree is not None:
        cmd += " -n " + newicktree

    if collapse:
        cmd += " --collapse"

    cmd += " -o " + str(outdir)

    subprocess.check_call(cmd, shell=True)

    typer.secho(f"Finished", fg=typer.colors.GREEN)


@tools_app.callback()
def items():
    """
    Contains various tools.
    """


def version_callback(value: bool):
    if value:
        typer.echo(f"SCGS CLI Version: {__version__}")
        raise typer.Exit()


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    verbose: bool = typer.Option(False, "--verbose/--silent", "-v/-s"),
    config: Path = typer.Option(
        "scgs.ini",
        exists=False,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        show_default=True,
    ),
    version: bool = typer.Option(None, "--version", callback=version_callback, is_eager=True),
):
    """
    Welcome to use gongyh/scgs pipeline!
    """
    # app_dir = typer.get_app_dir(APP_NAME)
    # typer.echo(f"This script runs at {app_dir}")
    typer.echo("Welcome to use gongyh/scgs pipeline!")
    if verbose:
        typer.echo("Will write verbose output")
        state["verbose"] = True

    if ctx.invoked_subcommand is None and config.exists():
        text = config.read_text()
        typer.echo("Check config file.")
        # typer.echo(f"Config file contents: \n{text}")


if __name__ == "__main__":
    app()
