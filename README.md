# gongyh/nf-core-scgs

**Single Cell Genome Sequencing data analysis pipeline**.

[![CircleCI](https://dl.circleci.com/status-badge/img/gh/gongyh/nf-core-scgs/tree/master.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/gongyh/nf-core-scgs/tree/master)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with podman](https://img.shields.io/badge/run%20with-podman-0dffed?labelColor=000000)](https://podman.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-ffb7ed?labelColor=000000&logo=anaconda)](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)

## Introduction

The pipeline is used for single cell genome sequencing data analysis and built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.

![Schematic image of scgs pipeline](scgs_pipeline.png)

## Quick start

Prerequisites: Git, Java 11 or later, Docker

```bash
## Install Nextflow
$ curl -s https://get.nextflow.io | bash
## Get the pipeline
$ git clone -b v2.0.1 https://github.com/gongyh/nf-core-scgs.git
## Test
$ ./nextflow run nf-core-scgs -profile test_local,docker
or $ ./nextflow run nf-core-scgs -profile test_local,podman
or $ APPTAINER_DISABLE_CACHE=true ./nextflow run nf-core-scgs -profile test_local,apptainer
or $ SINGULARITY_DISABLE_CACHE=true ./nextflow run nf-core-scgs -profile test_local,singularity
or $ ./nextflow run nf-core-scgs -profile test_local,conda
# for conda, add `disable_lockfile: true` to ~/.condarc or ~/.mambarc
```

## Documentation

The gongyh/nf-core-scgs pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
   - [Local installation](docs/configuration/local.md)
   - [Adding your own system](docs/configuration/adding_your_own.md)
   - [Reference genomes](docs/configuration/reference_genomes.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

## Related publications

[1] Jing, X., Gong, Y., Pan, H. et al. Single-cell Raman-activated sorting and cultivation (scRACS-Culture) for assessing and mining in situ phosphate-solubilizing microbes from nature. _ISME COMMUN_. 2, 106 (2022). [https://doi.org/10.1038/s43705-022-00188-3](https://doi.org/10.1038/s43705-022-00188-3)

[2] Jing X, Gong Y, Xu T, Meng Y, Han X, Su X, Wang J, Ji Y, Li Y, Jia Z, Ma B, Xu J. One-Cell Metabolic Phenotyping and Sequencing of Soil Microbiome by Raman-Activated Gravity-Driven Encapsulation (RAGE). _mSystems_. 2021 Jun 29;6(3):e0018121. doi: 10.1128/mSystems.00181-21. Epub 2021 May 27. PMID: [34042466](https://pubmed.ncbi.nlm.nih.gov/34042466/); PMCID: [PMC8269212](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc8269212/)

[3] Xu, T., Gong, Y., Su, X., Zhu, P., Dai, J., Xu, J., Ma, B., Phenome-Genome Profiling of Single Bacterial Cell by Raman-Activated Gravity-Driven Encapsulation and Sequencing. _Small_ 2020, 2001172. [https://doi.org/10.1002/smll.202001172](https://doi.org/10.1002/smll.202001172)

[4] Su, X., Gong, Y., Gou, H., Jing, X., Xu, T., Zheng, X., Chen, R., Li, Y., Ji, Y., Ma, B., Xu, J., Rational Optimization of Raman-Activated Cell Ejection and Sequencing for Bacteria. _Analytical Chemistry_, 2020. [https://doi.org/10.1021/acs.analchem.9b05345](https://doi.org/10.1021/acs.analchem.9b05345)

## Contact

gongyh/nf-core-scgs is developed by [Yanhai Gong](mailto:gongyh@qibebt.ac.cn), maintained by [Yanhai Gong](mailto:gongyh@qibebt.ac.cn) and [Shiqi Zhou](mailto:zhousq@qibebt.ac.cn). We look forward to receive your feedback, bug reports, or suggestions for the further development of this pipeline.

## [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This pipeline is open source under the MIT license, and integrates wonderful third-party softwares, which remain owned and copyrighted by their respective developers. Authors cannot be held legally or morally responsible for any consequences that may arise from using or misusing it.
