import os
import multiprocessing
from MGPGtools.utils.common import run


# Generate new gff and fasta files based on the gff file and core gene list
def gff2fasta(gff, refGenome, coreGenes, outdir):
    coregff = os.path.join(outdir, "core.gff")
    corefasta = os.path.join(outdir, "core.fasta")
    with open(gff, "r") as fgff, open(coregff, "w") as fcoregff:
        lines = fgff.read().strip().split("\n")
        for line in lines:
            if line[0] == "#":
                fcoregff.write(line + "\n")
                continue
            row = line.strip().split("\t")
            if row[2] == "CDS":
                gene = row[8].split(";")[0].replace("ID=", "")
                if gene in coreGenes:
                    fcoregff.write(line + "\n")
    extractCoreGenesCmd = [
        "gffread",
        coregff,
        "-g",
        refGenome,
        "-x",
        corefasta,
    ]
    isSuccess, stdout, stderr = run(extractCoreGenesCmd)


# Filter the nucmer results, extract gene IDs with a similarity greater than 80%
def filtNucmerResult(fasta, queryFa, geneLength, outdir):
    geneMapDict = {}
    nucmerCmd = [
        "nucmer",
        fasta,
        queryFa,
        "--delta",
        os.path.join(outdir, "out.delta"),
    ]
    run(nucmerCmd)
    showCoodsCmd = ["show-coords", os.path.join(outdir, "out.delta")]
    if_success, stdout, stderr = run(showCoodsCmd)
    i = 0
    for line in stdout.strip().split("\n"):
        i += 1
        if i < 6:
            continue
        row = line.strip().split()
        geneName = row[12].replace("_gene", "")
        geneMapDict[geneName] = {}
        geneMapDict[geneName]["ref_start"] = int(row[0])
        geneMapDict[geneName]["ref_end"] = int(row[1])
        if int(row[3]) > int(row[4]):
            geneMapDict[geneName]["query_start"] = int(row[4])
            geneMapDict[geneName]["query_end"] = int(row[3])
            geneMapDict[geneName]["orientation"] = "reverse"
        else:
            geneMapDict[geneName]["query_start"] = int(row[3])
            geneMapDict[geneName]["query_end"] = int(row[3])
            geneMapDict[geneName]["orientation"] = "forward"
        query_length = int(row[7])
        idy = float(row[9])
        if query_length * idy / (100 * geneLength[geneName]) > 0.8:
            geneMapDict[geneName]["assemblCore"] = True
        else:
            geneMapDict[geneName]["assemblCore"] = False
        geneMapDict[geneName]["refNode"] = row[11]
    return geneMapDict
