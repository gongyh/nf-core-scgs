# gongyh/nf-core-scgs: Databases

If running the pipeline in a local environment, you can download the databases of various tools using `bin/setup.sh` script.

```bash
docker run -v /data:/data -w /data -u $(id -u):$(id -g) -it gongyh/scgs /bin/bash
(scgs_py36) root@a9cea5b5a003:/data# /opt/nf-core-scgs/bin/setup.sh
```

## Supported databases

The following databases are supported by gongyh/scgs pipeline:

1. NCBI nt database (latest version)
2. Uniprot proteomes database (Reference_Proteomes_2017_07)
3. Kraken database (minikraken_20171101_8GB_dustmasked)
4. CheckM database (checkm_data_2015_01_16, no longer needed after v1.0-rc)
5. EggNOG v5 database (v5)
6. KOfam HMM database (latest)
7. Funannotate database (latest, only needed for eukaryotic microbes, e.g. fungi)
8. PlasmidFinder database (latest)
9. VirulenceFinder database (latest)
10. Pfam-A database (Pfam31.0)
11. EukCC database (V1.1, 20191023_1)
