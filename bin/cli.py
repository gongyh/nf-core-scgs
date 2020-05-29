#!/usr/bin/env python3.6

import time
from enum import Enum
from pathlib import Path
import typer
from Bio import SeqIO
import subprocess
import numpy as np

def real_split(fa, ann, level, out_dir):
    """
    Split by annotation
    """
    record_dict = SeqIO.to_dict(SeqIO.parse(fa, "fasta"))

    annCol = -1 # which colum corresponds to the specified level
    with open(ann) as fh:
        for line in fh:
            if line[0:2]=="##": # comment line, skip
                continue
            elif line[0:6]=="# name": # header line
                cl = line.strip().split("\t")
                for item in cl:
                    il = item.split(".")
                    if len(il)==3 and il[0]==level and il[1]=="t": # eg. order.t.12
                        annCol = int(il[2])-1
                continue
            cl = line.strip().split("\t")
            contig_id = cl[0]
            assert(annCol>0)
            annotation = cl[annCol]
            fname = annotation.replace(" ","_")+".fasta"
            with out_dir.joinpath(fname).open('a') as f:
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

@tools_app.command("scrs_pre")
def tools_scrs_preprocess(raw_dir: Path = typer.Option(
            ...,
            exists=True,
            file_okay=False,
            dir_okay=True,
            writable=False,
            readable=True,
            resolve_path=True,
           ),
           out_dir: Path = typer.Option(
            "./pre",
            exists=False,
            file_okay=False,
            dir_okay=True,
            writable=True,
            readable=True,
            resolve_path=True,
            show_default=True
           ),
           meta_table: Path = typer.Option(
            "./meta.txt",
            exists=True,
            file_okay=True,
            dir_okay=False,
            writable=False,
            readable=True,
            resolve_path=True,
            show_default=True
           ),
           out_prefix: str = typer.Option("SCRS", show_default=True)
      ):
    # first check input SCRS
    typer.echo(f"Checking input SCRS.")
    if raw_dir.exists() and raw_dir.is_dir():
        pass
    else:
        typer.echo(f"Please confirm {raw_dir} is exist and a directory.")
        typer.Exit()
    num_scrs = len(sorted(raw_dir.glob('*.txt')))
    if num_scrs == 0: # at least one spectrum needed
        typer.echo(f"No SCRS found, please check {raw_dir} .")
        typer.Exit()

    # then check meta table
    meta = np.loadtxt(meta_table, dtype=object, delimiter='\t')
    if len(meta.shape) != 2: # wrong shape
        typer.echo(f"Wrong format, please check {meta_table} .")
        typer.Exit()
    elif meta.shape[0] > num_scrs: #error
        typer.echo(f"Samples error, please check {meta_table} .")
        typer.Exit()
    elif meta.shape[0]>0 and meta.shape[1]<5: # not enough cols
        typer.echo(f"Not enough columns, please check {meta_table} .")
        typer.Exit()

    # call R script to process
    ret = subprocess.check_call(" ".join(['Rscript', str(Path(__file__).resolve().parent.joinpath('SCRS_preprocess.R')),
                                str(raw_dir), str(out_dir), str(meta_table), out_prefix]), shell=True)

    return ret

@tools_app.command("split")
def tools_split(results_dir: Path = typer.Option(
            "./results",
            exists=True,
            file_okay=False,
            dir_okay=True,
            writable=False,
            readable=True,
            resolve_path=True,
            show_default=True
           ), 
           level: TaxaLevels = typer.Option(TaxaLevels.genus, case_sensitive=False, show_default=True),
           output_dir: Path = typer.Option(
            "./split",
            exists=False,
            file_okay=False,
            dir_okay=True,
            writable=True,
            readable=True,
            resolve_path=True,
            show_default=True
           )
        ):
    typer.echo(f"Checking the existence of blob and spades results.")
    samples = []
    blob_dir = results_dir.joinpath('blob')
    if not blob_dir.is_dir():
        typer.echo(f"Taxa annotations not found, please check {blob_dir} .")
        typer.Exit()
    else: # check subdir
        for child in blob_dir.iterdir():
            if child.is_dir(): # sample
                samples.append(child.parts[-1])
    typer.echo(f"INFO: {len(samples)} samples detected.")
    spades_dir = results_dir.joinpath('spades')
    if not spades_dir.is_dir():
        typer.echo(f"Spades assemblies not found, please check {spades_dir} .")
        typer.Exit()
    else: # check the existence of all genome assemblies
        for sample in samples:
            ass = spades_dir.joinpath(sample+'.ctg200.fasta')
            if not ass.exists():
                typer.echo(f"genome assembly file {ass} not found.")
                typer.Exit()
    typer.echo(f"Done!")

    # try to create outdir
    try:
        output_dir.mkdir(exist_ok=True)
    except:
        typer.echo(f"Can not create output directory.")
        typer.Exit()

    # process split one by one
    with typer.progressbar(samples, label="Processing") as progress:
        for sample in progress:
            #typer.echo(sample)
            fa = spades_dir.joinpath(sample+'.ctg200.fasta')
            out_subdir = output_dir.joinpath(sample+'_'+level)
            out_subdir.mkdir(exist_ok=True) # create subdir to store fastas
            blob_sub = blob_dir.joinpath(sample)
            blob_table = None
            for child in blob_sub.iterdir():
                if child.match(sample+'.blobDB*.table.txt'):
                    blob_table = child
            if blob_table is None:
                typer.echo(f"Can not find annotation table for sample {sample}.")
                typer.Exit()
            # perform split
            real_split(fa, blob_table, level, out_subdir)
        typer.echo(f"\nFinished.")

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
def main(ctx: typer.Context, verbose: bool = typer.Option(False, "--verbose/--silent", "-v/-s"),
         config: Path = typer.Option(
           "scgs.ini",
           exists=False,
           file_okay=True,
           dir_okay=False,
           writable=False,
           readable=True,
           resolve_path=True,
           show_default=True
           ),
         version: bool = typer.Option(None, "--version", callback=version_callback, is_eager=True),
        ):
    """
    Welcome to use gongyh/scgs pipeline!
    """
    app_dir = typer.get_app_dir(APP_NAME)
    typer.echo(f"This script runs at {app_dir}")
    if verbose:
        typer.echo("Will write verbose output")
        state["verbose"] = True

    if ctx.invoked_subcommand is None and config.exists():
        text = config.read_text()
        typer.echo("Check config file.")
        typer.echo(f"Config file contents: \n{text}")


if __name__ == "__main__":
    app()

