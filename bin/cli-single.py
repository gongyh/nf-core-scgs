#!/usr/bin/env python3

import os
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
            bacAnnotation = cl[annCol].replace("/", "_")
            eukAnnotation = cl[eukAnnCol].replace("/", "_")
            if "_" in bacAnnotation:
                bacAnnotation = bacAnnotation.split("_")[0]
            if "_" in eukAnnotation:
                eukAnnotation = eukAnnotation.split("_")[0]
            superkingdom = cl[eukCol]
            bac_gname = bacAnnotation.replace(" ", "_") + ".gids"
            bac_kname = bacAnnotation.replace(" ", "_") + ".ko"
            euk_gname = eukAnnotation.replace(" ", "_") + ".gids"
            euk_kname = eukAnnotation.replace(" ", "_") + ".ko"
            euk_kofh = euk_out_dir.joinpath(euk_kname).open("a")
            euk_gname_path = euk_out_dir.joinpath(euk_gname)
            bac_kofh = bac_out_dir.joinpath(bac_kname).open("a")
            bac_gname_path = bac_out_dir.joinpath(bac_gname)
            with bac_gname_path.open("a") as f1, euk_gname_path.open("a") as f2:
                if superkingdom == "Eukaryota":
                    if ctg_id_short in ctg_genes.keys():
                        for g in ctg_genes[ctg_id_short]:
                            f2.write(g + "\n")
                            if g in gene_ko.keys():
                                for v in gene_ko[g]:
                                    euk_kofh.write(g + "\t" + v + "\n")
                else:
                    if ctg_id_short in ctg_genes.keys():
                        for g in ctg_genes[ctg_id_short]:
                            f1.write(g + "\n")
                            if g in gene_ko.keys():
                                for v in gene_ko[g]:
                                    bac_kofh.write(g + "\t" + v + "\n")
            bac_kofh.close()
            euk_kofh.close()
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
    sample: str = typer.Option(help="Samples Name"),
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


    blob_dir = results_dir.joinpath("blob")
    if not blob_dir.is_dir():
        typer.secho(
            f"Taxa annotations not found, please check {blob_dir} .",
            fg=typer.colors.RED,
        )
        raise typer.Abort()
    if sample not in blob_dir.listdir():
        typer.secho(
            f"Sample taxa annotations not found, please check {blob_dir} .",
            fg=typer.colors.RED,
        )
        raise typer.Abort()
    typer.echo(f"INFO: {sample} detected.")

    spades_dir = results_dir.joinpath("spades")
    if not spades_dir.is_dir():
        typer.secho(
            f"Spades assemblies not found, please check {spades_dir} .",
            fg=typer.colors.RED,
        )
        raise typer.Abort()
    else:  # check the existence of all genome assemblies
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
