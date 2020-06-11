#!/usr/bin/env python3.6

import time
from enum import Enum
from pathlib import Path
import typer
from Bio import SeqIO
import subprocess
import numpy as np

def real_split(fa, ann, level, out_dir, gff=None, ko_file=None):
    """
    Split by annotation
    """
    record_dict = SeqIO.to_dict(SeqIO.parse(fa, "fasta"))

    ctg_genes = {}
    if gff is None: # skip splitting gene ids
        pass
    else: # parse the gff file
        with open(gff) as gfh:
            for line in gfh:
                cl = line.strip().split("\t")
                if len(cl) >= 9: # gene annotation line
                    ctg_id = cl[0]
                    anno8 = cl[8].split(";")[0]
                    if anno8.startswith("ID=") and anno8.endswith("_gene"): # correct line
                        gid = anno8[3:len(anno8)-5]
                        if ctg_id in ctg_genes.keys(): # add the new id
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
                if len(cl) == 2: # correct annotation
                    if cl[0] not in gene_ko.keys():
                        gene_ko[cl[0]] = [cl[1]]
                    else:
                        gene_ko[cl[0]].append(cl[1])

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
            ctg_id_short = contig_id.split("_length_")[0]
            annotation = cl[annCol]
            gname = annotation.replace(" ","_")+".gids"
            kname = annotation.replace(" ","_")+".ko"
            kofh = out_dir.joinpath(kname).open('a')
            with out_dir.joinpath(gname).open('a') as f:
                if ctg_id_short in ctg_genes.keys():
                    for g in ctg_genes[ctg_id_short]:
                        f.write(g+'\n')
                        if g in gene_ko.keys():
                            for v in gene_ko[g]:
                                kofh.write(g+"\t"+v+"\n")
            kofh.close()
            assert(annCol>0)
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

@tools_app.command("scrs_filter", help="Remove low quality SCRS")
def tools_scrs_filter(raw_dir: Path = typer.Option(
            ...,
            exists=True,
            file_okay=False,
            dir_okay=True,
            writable=False,
            readable=True,
            resolve_path=True,
           ),
           out_dir: Path = typer.Option(
            "./good/",
            exists=False,
            file_okay=False,
            dir_okay=True,
            writable=True,
            readable=True,
            resolve_path=True,
            show_default=True
           )
         ):
    typer.echo(f"Checking input SCRS.")
    if raw_dir.exists() and raw_dir.is_dir():
        pass
    else:
        typer.secho(f"Please confirm {raw_dir} is exist and a directory.", fg=typer.colors.RED)
        raise typer.Abort()
    num_scrs = len(sorted(raw_dir.glob('*.txt')))
    if num_scrs == 0: # at least one spectrum needed
        typer.secho(f"No SCRS found, please check {raw_dir} .", fg=typer.colors.RED)
        raise typer.Abort()

    # call R script to process
    subprocess.check_call(" ".join(['Rscript', str(Path(__file__).resolve().parent.joinpath('SCRS_filter.R')),
                          str(raw_dir), str(out_dir)]), shell=True)

    typer.secho(f"Finished", fg=typer.colors.GREEN)

@tools_app.command("scrs_pre", help="Preprocess for SCRS")
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
            "./pre/",
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
        typer.secho(f"Please confirm {raw_dir} is exist and a directory.", fg=typer.colors.RED)
        raise typer.Abort()
    num_scrs = len(sorted(raw_dir.glob('*.txt')))
    if num_scrs == 0: # at least one spectrum needed
        typer.secho(f"No SCRS found, please check {raw_dir} .", fg=typer.colors.RED)
        raise typer.Abort()

    # then check meta table
    meta = np.loadtxt(meta_table, dtype=object, delimiter='\t')
    #typer.echo(f"Size of meta table: {meta.shape}")
    if len(meta.shape) != 2: # wrong shape
        typer.secho(f"Wrong format, please check {meta_table} .", fg=typer.colors.RED)
        raise typer.Abort()
    elif meta.shape[0] <= 1: # not choose any SCRS
        typer.secho(f"Samples error, please check {meta_table} .", fg=typer.colors.RED)
        raise typer.Abort()
    elif meta.shape[1]<5: # not enough cols
        typer.secho(f"Not enough columns, please check {meta_table} .", fg=typer.colors.RED)
        raise typer.Abort()

    # call R script to process
    subprocess.check_call(" ".join(['Rscript', str(Path(__file__).resolve().parent.joinpath('SCRS_preprocess.R')),
                          str(raw_dir), str(out_dir), str(meta_table), out_prefix]), shell=True)

    typer.secho(f"Finished", fg=typer.colors.GREEN)

@tools_app.command("scrs_snr", help="Signal-Noise ratio calculation")
def tools_scrs_snr(scgs_csv: Path = typer.Option(
            ...,
            exists=True,
            file_okay=True,
            dir_okay=False,
            writable=False,
            readable=True,
            resolve_path=True,
            show_default=True
           ),
           out_dir: Path = typer.Option(
            "./stats/",
            exists=False,
            file_okay=False,
            dir_okay=True,
            writable=True,
            readable=True,
            resolve_path=True,
            show_default=True
           )
         ):
    # first check input SCRS
    typer.echo(f"Checking input SCRS.")
    if scgs_csv.exists() and scgs_csv.is_file():
        pass
    else:
        typer.secho(f"Please confirm {scgs_csv} is exist and a CSV file.", fg=typer.colors.RED)
        raise typer.Abort()

    # call R script to process
    subprocess.check_call(" ".join(['Rscript', str(Path(__file__).resolve().parent.joinpath('SCRS_snr.R')),
                          str(scgs_csv), str(out_dir)]), shell=True)

    typer.secho(f"Finished", fg=typer.colors.GREEN)

@tools_app.command("scrs_cdr", help="C-H / C-D peaks calculation")
def tools_scrs_cdr(scgs_csv: Path = typer.Option(
            ...,
            exists=True,
            file_okay=True,
            dir_okay=False,
            writable=False,
            readable=True,
            resolve_path=True,
            show_default=True
           ),
           out_dir: Path = typer.Option(
            "./stats/",
            exists=False,
            file_okay=False,
            dir_okay=True,
            writable=True,
            readable=True,
            resolve_path=True,
            show_default=True
           )
         ):
    # first check input SCRS
    typer.echo(f"Checking input SCRS.")
    if scgs_csv.exists() and scgs_csv.is_file():
        pass
    else:
        typer.secho(f"Please confirm {scgs_csv} is exist and a CSV file.", fg=typer.colors.RED)
        raise typer.Abort()

    # call R script to process
    subprocess.check_call(" ".join(['Rscript', str(Path(__file__).resolve().parent.joinpath('SCRS_cdr.R')),
                          str(scgs_csv), str(out_dir)]), shell=True)

    typer.secho(f"Finished", fg=typer.colors.GREEN)

@tools_app.command("scrs_rarefy", help="Sample depth for CDR")
def tools_scrs_rarefy(scgs_cdr: Path = typer.Option(
            ...,
            exists=True,
            file_okay=True,
            dir_okay=False,
            writable=False,
            readable=True,
            resolve_path=True,
            show_default=True
           ),
           out_dir: Path = typer.Option(
            "./rarefy/",
            exists=False,
            file_okay=False,
            dir_okay=True,
            writable=True,
            readable=True,
            resolve_path=True,
            show_default=True
           )
         ):
    # first check input SCRS
    typer.echo(f"Checking input SCRS.")
    if scgs_cdr.exists() and scgs_cdr.is_file():
        pass
    else:
        typer.secho(f"Please confirm {scgs_cdr} is exist and a CSV file.", fg=typer.colors.RED)
        raise typer.Abort()

    # call R script to process
    subprocess.check_call(" ".join(['Rscript', str(Path(__file__).resolve().parent.joinpath('SCRS_cdr_stats.R')),
                          str(scgs_cdr), str(out_dir)]), shell=True)
    stats_out = out_dir.joinpath('CD_ratio_cummean_accum_df_by_Time.txt')
    if stats_out.exists() and stats_out.is_file():
        subprocess.check_call(" ".join(['Rscript', str(Path(__file__).resolve().parent.joinpath('SCRS_cdr_rarefy.R')),
                          str(stats_out), str(out_dir)]), shell=True)
    else:
        typer.secho(f"Can not perform CDR statistics for {scgs_cdr}.", fg=typer.colors.RED)
        raise typer.Abort()

    typer.secho(f"Finished", fg=typer.colors.GREEN)


@tools_app.command("scrs_pipeline", help="Pick->normalize->S/N->CDR->rarefy")
def tools_scrs_pipeline(raw_dir: Path = typer.Option(
            ...,
            exists=True,
            file_okay=False,
            dir_okay=True,
            writable=False,
            readable=True,
            resolve_path=True,
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
           out_dir: Path = typer.Option(
            "./scrs/",
            exists=False,
            file_okay=False,
            dir_okay=True,
            writable=True,
            readable=True,
            resolve_path=True,
            show_default=True
           ),
           cdr: bool = typer.Option(False, "--cdr", "-c"),
           rarefy: bool = typer.Option(False, "--rarefy", "-r")
         ):
    out_dir.mkdir(exist_ok=True)
    typer.secho(f"INFO: filter low quality SCRS.", fg=typer.colors.GREEN)
    tools_scrs_filter(raw_dir, out_dir.joinpath('good'))
    typer.secho(f"INFO: preprocess SCRS.", fg=typer.colors.GREEN)
    tools_scrs_preprocess(out_dir.joinpath('good'), out_dir.joinpath('pre'), meta_table,'SCRS')
    typer.secho(f"INFO: Calc SNR for SCRS.", fg=typer.colors.GREEN)
    tools_scrs_snr(out_dir.joinpath('pre','Cells_bg_baseline_zero_scale.csv'), out_dir.joinpath('stats'))
    typer.secho(f"INFO: Calc CDR for SCRS.", fg=typer.colors.GREEN)
    if cdr:
        tools_scrs_cdr(out_dir.joinpath('pre','Cells_bg_baseline_zero_scale.csv'), out_dir.joinpath('stats'))
        if rarefy: # only when cdr ready
            typer.secho(f"INFO: Rarefy based on CDR.", fg=typer.colors.GREEN)
            tools_scrs_rarefy(out_dir.joinpath('stats','CDR.txt'), out_dir.joinpath('rarefy'))
    typer.secho(f"Done!", fg=typer.colors.GREEN)

####################################################################################

@tools_app.command("scgs_split", help="Split assembly genome according to taxa")
def tools_split(results_dir: Path = typer.Option(
            "./results/",
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
            "./split/",
            exists=False,
            file_okay=False,
            dir_okay=True,
            writable=True,
            readable=True,
            resolve_path=True,
            show_default=True
           ),
           force: bool = typer.Option(False, "--force", "-f")
        ):
    typer.echo(f"Checking the existence of blob, prokka and spades results.")
    samples = []

    blob_dir = results_dir.joinpath('blob')
    if not blob_dir.is_dir():
        typer.secho(f"Taxa annotations not found, please check {blob_dir} .", fg=typer.colors.RED)
        raise typer.Abort()
    else: # check subdir
        for child in blob_dir.iterdir():
            if child.is_dir(): # sample
                samples.append(child.parts[-1])
    typer.echo(f"INFO: {len(samples)} samples detected.")

    spades_dir = results_dir.joinpath('spades')
    if not spades_dir.is_dir():
        typer.secho(f"Spades assemblies not found, please check {spades_dir} .", fg=typer.colors.RED)
        raise typer.Abort()
    else: # check the existence of all genome assemblies
        for sample in samples:
            ass = spades_dir.joinpath(sample+'.ctg200.fasta')
            if not ass.exists():
                typer.secho(f"genome assembly file {ass} not found.", fg=typer.colors.RED)
                raise typer.Abort()

    prokka_dir = results_dir.joinpath('prokka')
    anno_exist = False
    if prokka_dir.exists() and prokka_dir.is_dir(): # exist the dir, ok
        typer.echo(f"INFO: gene annotations found, will also split gene ids.")
        anno_exist = True
    else:
        typer.echo(f"INFO: gene annotations not found, will skip splitting gene ids.")
        anno_exist = False

    kofam_dir = results_dir.joinpath('kofam')
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
            output_dir.rename(str(output_dir)+"%d"%(int(time.time())))
            output_dir.mkdir(exist_ok=False)
        else:
            output_dir.mkdir(exist_ok=True)
    except:
        typer.secho(f"Can not create output directory.", fg=typer.colors.RED)
        raise typer.Abort()

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
                typer.secho(f"Can not find annotation table for sample {sample}.", fg=typer.colors.RED)
                raise typer.Abort()
            gff = None
            sample_anno = prokka_dir.joinpath(sample)
            if anno_exist and sample_anno.exists() and sample_anno.is_dir():
                sample_gff = sample_anno.joinpath(sample+".gff")
                if sample_gff.exists() and sample_gff.is_file(): # find annotation gff file
                    gff = sample_gff
            ko_file = None
            if ko_exist:
                ko_file = kofam_dir.joinpath(sample+"_KOs_mapper.txt")
            # perform split
            real_split(fa, blob_table, level, out_subdir, gff, ko_file)
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
    #app_dir = typer.get_app_dir(APP_NAME)
    #typer.echo(f"This script runs at {app_dir}")
    typer.echo("Welcome to use gongyh/scgs pipeline!")
    if verbose:
        typer.echo("Will write verbose output")
        state["verbose"] = True

    if ctx.invoked_subcommand is None and config.exists():
        text = config.read_text()
        typer.echo("Check config file.")
        #typer.echo(f"Config file contents: \n{text}")


if __name__ == "__main__":
    app()

