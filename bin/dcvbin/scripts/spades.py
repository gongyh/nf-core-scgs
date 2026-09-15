import os
import subprocess
import argparse

def run_command(cmd):
    print(f"📦 Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def run_pipeline(input_folder, threads, contig_min_len=2000):
    files = os.listdir(input_folder)
    fastq1 = next((f for f in files if f.endswith('_1.fastq.gz')), None)
    fastq2 = next((f for f in files if f.endswith('_2.fastq.gz')), None)

    if not fastq1 or not fastq2:
        print("❌ Cannot find _1.fastq.gz and _2.fastq.gz files.")
        return

    base_name = os.path.basename(fastq1).split('_')[0]
    fq1_path = os.path.join(input_folder, fastq1)
    fq2_path = os.path.join(input_folder, fastq2)

    output_dir = os.path.join(input_folder, "assembly_out")
    contigs_path = os.path.join(output_dir, "contigs.fasta")
    filtered_contigs = os.path.join(input_folder, "contigs_over2k.fasta")

    print("🔧 Step 1: Running SPAdes...")
    run_command(
        f"spades.py --meta -1 {fq1_path} -2 {fq2_path} -o {output_dir} -t {threads}"
    )

    print("🔧 Step 2: Filtering contigs > {} bp...".format(contig_min_len))
    run_command(
        f"""awk '/^>/ {{if (seqlen>={contig_min_len}) print seqname"\\n"seq; seqname=$0; seq=""; seqlen=0; next}}
               {{seq=seq$0; seqlen+=length($0)}}
               END {{if (seqlen>={contig_min_len}) print seqname"\\n"seq}}' \
               {contigs_path} > {filtered_contigs}"""
    )

    print("🔧 Step 3: BWA index...")
    run_command(f"bwa index {filtered_contigs}")

    sam_out = os.path.join(input_folder, f"{base_name}.sam")
    bam_out = os.path.join(input_folder, f"{base_name}.bam")
    sorted_bam = os.path.join(input_folder, f"{base_name}_sorted.bam")

    print("🔧 Step 4: BWA mem...")
    run_command(
        f"bwa mem {filtered_contigs} {fq1_path} {fq2_path} -t {threads} -o {sam_out}"
    )

    print("🔧 Step 5: SAM to BAM...")
    run_command(
        f"samtools view -F 3584 -b --threads {threads} {sam_out} > {bam_out}"
    )

    print("🔧 Step 6: Sorting BAM...")
    run_command(
        f"samtools sort {bam_out} -o {sorted_bam} -@ {threads}"
    )

    print("✅ All steps completed! Output BAM: {}".format(sorted_bam))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="One-click SPAdes + BWA + BAM pipeline")
    parser.add_argument("-i", "--input", required=True, help="Input folder containing fastq.gz files")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads")
    args = parser.parse_args()

    run_pipeline(args.input, args.threads)
