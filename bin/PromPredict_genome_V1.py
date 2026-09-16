#!/usr/bin/env python

import argparse
import pexpect
import sys

def run_prom_predict_with_args():
    # ----------------------
    # Step 1: Parse command-line arguments with argparse
    # ----------------------
    parser = argparse.ArgumentParser(
        description="Wrapper script to run PromPredict_genome_V1 with interactive parameters via command line."
    )

    # Define the three required parameters (match PromPredict's interactive inputs)
    parser.add_argument(
        "--genome_fasta",
        type=str,
        required=True,
        help="Path to input genome FASTA file"
    )

    parser.add_argument(
        "--window_size",
        type=str,
        required=False,
        default="100",
        help="Window size for the upstream region to be compared. Advised size: 100nt (default: 100)"
    )

    parser.add_argument(
        "--gc",
        type=str,
        default="default",
        required=False,
        help="GC-content of the whole genome. Use the 1000nt fragement %%GC to set this cut-off value, type 'default' (default: default)"
    )

    # Parse arguments
    args = parser.parse_args()


    # ----------------------
    # Step 2: Use pexpect to run PromPredict_genome_V1 interactively
    # ----------------------
    # Replace with the actual path to PromPredict_genome_V1 (e.g., "/opt/tools/PromPredict_genome_V1")
    PROM_PREDICT_PATH = "PromPredict_genome_V1"

    try:
        # Start the PromPredict program
        print(f"Starting PromPredict_genome_V1...")
        child = pexpect.spawn(
            PROM_PREDICT_PATH,
            timeout=60  # Timeout (seconds) for each interaction step
        )

        # ----------------------
        # Send parameters to interactive prompts
        # ----------------------

        # 1st parameter: Genome FASTA path
        child.expect("Enter the Input genome File Name:", timeout=30)
        child.sendline(args.genome_fasta)
        print(f"Sent genome FASTA path: {args.genome_fasta}")

        # 2nd parameter: E1 region window size
        child.expect("Enter the E1 region window size", timeout=30)
        child.sendline(args.window_size)
        print(f"Sent E1 window size: {args.window_size}")

        # 3rd parameter: GC%
        child.expect("Enter the whole genome GC content:", timeout=30)
        child.sendline(args.gc)
        print(f"Sent GC content: {args.gc}")

        # ----------------------
        # Wait for program completion
        # ----------------------
        print("PromPredict_genome_V1 is running...")
        child.expect(pexpect.EOF, timeout=None)  # Wait indefinitely for program to finish
        print("\nPromPredict_genome_V1 completed successfully!")

    except pexpect.TIMEOUT:
        print("\nError: Timeout waiting for program prompt. Check if prompts match PromPredict's output!")
    except pexpect.EOF:
        print("\nError: PromPredict_genome_V1 exited unexpectedly.")
    finally:
        child.close()  # Ensure the child process is terminated


if __name__ == "__main__":
    run_prom_predict_with_args()
