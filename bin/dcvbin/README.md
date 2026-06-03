# DCVBin: A Novel Metagenomic Binning Tool

## Overview

**DCVBin** is an innovative metagenomic binning tool that combines **compositional features** and **semantic features** to improve the recovery of high-quality Metagenome-Assembled Genomes (MAGs). By leveraging a DNA language model (DNABERT_S) to extract semantic features and integrating them with compositional data, DCVBin addresses binning challenges, especially in cases with limited coverage information. This tool enhances binning performance and provides better recovery of high-quality MAGs from metagenomic data.

The model provides two optional analysis modes:

### Standard Mode
Uses the pre-trained DNABERT model to extract sequence features, combined with k-mer frequency and sample-specific coverage features. After feature fusion via Variational Autoencoders (VAE), two-step adaptive clustering is performed to bin the sequences.

### Enhanced Mode
Builds on the standard process by adding domain-adaptive pretraining of DNABERT, fine-tuning it on your data. This improves feature extraction for specific microbiomes, yielding better results in identifying new species, but requires more computational resources.

Both modes automatically complete the entire pipeline, from raw data to quality assessment. The enhanced mode incorporates domain adaptation training, significantly improving the identification of new species, but requiring more computational resources. Users can choose the most suitable analysis mode based on their specific research goals and available computational resources.

## Installation

To install **DCVBin**, follow these steps:

1. **Clone the repository** to your local machine:

    ```bash
    git clone https://github.com/dengdengf/dcvbin.git
    cd DCVBin
    ```

2. **Create and activate a conda environment** using the provided `.yaml` file in the `envs` directory:

    ```bash
    conda env create -f envs/checkm2.yaml
    conda env create -f envs/copygen.yaml
    conda env create -f envs/dnaberts.yaml
    conda create --name vambnew --file envs/vambnew.lock
    conda activate dcvbin_env
    ```


## Usage

### 1. **Basic Mode (without CPT)**

In the **basic mode**, you will run the `my_method.sh` script to process your input FASTA files and generate MAGs.

- **Prepare your data**: Ensure that your input data is in the `.fasta` format.Ensure your input data is in FASTA format. If the data is in a different format, consider converting it to FASTA first.
  
- **Run the script**:

    ```bash
    bash my_method.sh input_dir output_dir
    ```

    - `input_dir`: Directory containing your `.fasta` files.
    - `output_dir`: Directory where the output MAGs will be saved.

### 2. **CPT Mode (with Continued Pretraining)**

For the **CPT mode**, use the `my_method_fur.sh` script to run the process with **Continued Pretraining (CPT)** of the DNA language model.

- **Run the script**:

    ```bash
    bash my_method_fur.sh input_dir output_dir
    ```

    - `input_dir`: Directory containing your `.fasta` files.
    - `output_dir`: Directory where the output MAGs will be saved.

The `DNABERT_S-main` directory is too large for GitHub and has been uploaded to Figshare:

**Download Link**: (https://figshare.com/articles/figure/DCVBin_-_Training_Data_and_Models_for_Metagenomic_Binning/30400258)

If you need to use CPT mode, you must:

1. **Download** the compressed file from the above link
2. **Extract** the archive to your project root directory:
   ```bash
   tar -xzf DNABERT_S-main.tar.gz

### Example:

To run in **basic mode**:

```bash
bash my_method.sh ./data ./data/output
```

To run in CPT mode:

```bash 
bash my_method_fur.sh ./data ./data/output
```

Contact

For any questions, concerns, or issues with DCVBin, please contact LiuYF via email at 17837410609@163.com, or create an issue on the [GitHub Issues page](https://github.com/dengdengf/dcvbin/issues).

Thank you for using DCVBin!
