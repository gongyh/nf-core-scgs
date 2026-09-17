# Database Configuration

`gongyh/nf-core-scgs` can download the databases used by its SCGS and
MINIMETA workflows. Run the dedicated database workflow once in a persistent
location, then reuse the downloaded directories across analyses.

Database downloads require network access and can consume substantial storage.
Run the preparation workflow on storage that is accessible to the compute
environment used for the analysis workflow.

## Prepare Databases

Use `--prepare_databases` to select the database workflow:

```bash
nextflow run gongyh/nf-core-scgs \
    --prepare_databases \
    --db_type checkm2,mmseqs,kofam,eggnog \
    -profile docker
```

The workflow writes prepared resources to `./databases` by default and records
the software versions used in `pipeline_info/`.

To prepare every supported database, omit `--db_type` or set it to `all`:

```bash
nextflow run gongyh/nf-core-scgs \
    --prepare_databases \
    --db_type all \
    -profile docker
```

`--db_type` accepts a comma-separated, case-insensitive list. Valid values
are:

```text
mmseqs, checkm2, kofam, eggnog, kraken2, gtdb, blob, metabuli, genomad, nt, all
```

Do not use this workflow for every analysis by default. Prepare only the
resources needed for the options you intend to enable.

## Database Map

The table below maps each preparation selector to its output and the pipeline
parameter that consumes it. Paths are shown relative to the database output
directory.

| `--db_type` value | Prepared resource                           | Use it with                                                         | Workflow          |
| ----------------- | ------------------------------------------- | ------------------------------------------------------------------- | ----------------- |
| `checkm2`         | `checkm2_db/`                               | `--checkm2_db /path/to/checkm2_db`                                  | SCGS and MINIMETA |
| `mmseqs`          | `mmseqs_db/`                                | `--mmseqs_db /path/to/mmseqs_db`                                    | MINIMETA          |
| `kofam`           | `kofam_db/profiles/` and `kofam_db/ko_list` | `--kofam_profile /path/to/profiles --kofam_kolist /path/to/ko_list` | SCGS and MINIMETA |
| `eggnog`          | `eggnog_db/`                                | `--eggnog_db /path/to/eggnog_db`                                    | SCGS and MINIMETA |
| `kraken2`         | `kraken2_db/`                               | `--kraken2_db /path/to/kraken2_db`                                  | SCGS              |
| `gtdb`            | `gtdb_db/`                                  | `--gtdb /path/to/gtdb_db`                                           | SCGS              |
| `blob`            | `blob_db/nodesDB.txt`                       | `--blob_db /path/to/blob_db/nodesDB.txt`                            | SCGS              |
| `metabuli`        | `metabuli_db/`                              | `--metabuli_db /path/to/metabuli_db`                                | MINIMETA          |
| `genomad`         | `db/`                                       | `--genomad_db /path/to/db`                                          | SCGS              |
| `nt`              | `nt_db/`                                    | `--nt_db /path/to/nt_db`                                            | SCGS              |

`checkm2` enables CheckM2 quality assessment. In MINIMETA,
`--run_cooccurrence_checkm` additionally runs CheckM2 on co-occurrence bins.

`mmseqs` provides contig taxonomy for MINIMETA and can improve SemiBin2
binning. `metabuli` provides taxonomy for the optional TaxVAMB branch.

KOfam annotation requires both `--kofam_profile` and `--kofam_kolist`.
EggNOG annotation requires `--eggnog_db`. In both cases, the corresponding
workflow flag (`--kofam` or `--eggnog`) must remain enabled.

## Reuse Databases With A Config File

For repeatable runs, place absolute paths in a project or site config file:

```nextflow
params {
    checkm2_db    = '/shared/scgs-databases/checkm2_db'
    mmseqs_db     = '/shared/scgs-databases/mmseqs_db'
    metabuli_db   = '/shared/scgs-databases/metabuli_db'
    kraken2_db    = '/shared/scgs-databases/kraken2_db'
    gtdb          = '/shared/scgs-databases/gtdb_db'
    blob_db       = '/shared/scgs-databases/blob_db/nodesDB.txt'
    genomad_db    = '/shared/scgs-databases/db'
    nt_db         = '/shared/scgs-databases/nt_db'
    kofam_profile = '/shared/scgs-databases/kofam_db/profiles'
    kofam_kolist  = '/shared/scgs-databases/kofam_db/ko_list'
    eggnog_db     = '/shared/scgs-databases/eggnog_db'
}
```

Load the configuration when starting an analysis:

```bash
nextflow run gongyh/nf-core-scgs \
    --minimeta \
    --reads 'reads/*_R{1,2}.fastq.gz' \
    -c databases.config \
    -profile docker
```

The same configuration can be used with SCGS runs. A resource is used only
when the related analysis option is enabled; unused database paths do not
trigger their associated analysis steps.

## Resources Not Prepared Here

The preparation workflow does not download every optional SCGS resource. You
must obtain and configure the following independently when needed:

| Resource                                      | Parameter                          | Used for                           |
| --------------------------------------------- | ---------------------------------- | ---------------------------------- |
| Kraken1 database                              | `--kraken1_db`                     | ACDC                               |
| Krona taxonomy file                           | `--krona_db`                       | Offline Krona reports              |
| UniProt protein database and taxonomy mapping | `--uniprot_db`, `--uniprot_taxids` | DIAMOND and BlobTools annotation   |
| Trusted Prokka proteins                       | `--prokka_proteins`                | Prokka annotation                  |
| Bakta database                                | `--bakta_db`                       | Bakta annotation                   |
| EukCC database                                | `--eukcc_db`                       | Eukaryotic completeness assessment |
| MGPG database                                 | `--mgpg_db`                        | Pangenome analysis                 |
| DNABERT-S model directory                     | `--DNABERTS_dir`                   | MINIMETA DCVBIN integration        |

See the [SCGS usage guide](../usage.md) and the
[MINIMETA workflow guide](../../MINIMETA.md) for the analysis options that
activate these resources.

## Operational Notes

- Keep databases outside Nextflow `work/` directories. The workflow output is
  the reusable copy; task work directories can be cleaned safely after a
  successful run.
- Use absolute paths in shared configuration files so that compute nodes and
  container runtimes resolve the same resources.
- With Docker, Podman, Singularity, or Apptainer, ensure the database parent
  directory is mounted or otherwise visible to task containers.
- Re-run a specific `--db_type` value when you intentionally need to refresh
  that resource.
