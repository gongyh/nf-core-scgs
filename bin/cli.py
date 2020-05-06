#!/usr/bin/env python3.6

import time
from enum import Enum
from pathlib import Path
import typer

class NeuralNetwork(str, Enum):
    simple = "simple"
    conv = "conv"
    lstm = "lstm"
    

state = {"verbose": False}
__version__ = "0.1.0"
APP_NAME = "my-super-cli-app"

app = typer.Typer()
items_app = typer.Typer()
app.add_typer(items_app, name="items", help="Manage items.")

@items_app.command("create")
def items_create(item: str):
    typer.echo(f"Creating item: {item}")

@items_app.callback()
def items():
    """
    Manage items in the app.
    """
    users = ["Camila", "Rick", "Morty"]

    with typer.progressbar(users, label="Processing") as progress:
        for user in progress:
            time.sleep(1)
        typer.echo(f"\nFinished.")


def version_callback(value: bool):
    if value:
        typer.echo(f"Awesome CLI Version: {__version__}")
        raise typer.Exit()

@app.command("say_hello", help="Say hello to NAME.")
def hello(name: str = typer.Option(..., "--name", "-n"),
          lastname: str = typer.Option("", help="this option does this and that."),
          version: bool = typer.Option(None, "--version", callback=version_callback, is_eager=True),
         ):
    """
    Say hello to NAME.
    """
    typer.echo(f"Hello {name}")


@app.command()
def goodbye(name: str, formal: bool = False):
    if formal:
        typer.echo(f"Goodbye Ms. {name}. Have a good day.")
    else:
        typer.echo(f"Bye {name}!")


@app.callback(invoke_without_command=True)
def main(ctx: typer.Context, verbose: bool = typer.Option(False, "--verbose/--silent", "-v/-s"),
         network: NeuralNetwork = typer.Option(NeuralNetwork.simple, case_sensitive=False, show_default=True),
         config: Path = typer.Option(
           "~/.scgs.ini",
           exists=True,
           file_okay=True,
           dir_okay=False,
           writable=False,
           readable=True,
           resolve_path=True,
           show_default=True
           )
        ):
    """
    Welcome to use gongyh/scgs pipeline!
    """
    app_dir = typer.get_app_dir(APP_NAME)
    if verbose:
        typer.echo("Will write verbose output")
        state["verbose"] = True

    if ctx.invoked_subcommand is None:
        text = config.read_text()
        typer.echo("Check config file.")
        typer.echo(f"Config file contents: \n{text}")


if __name__ == "__main__":
    app()

