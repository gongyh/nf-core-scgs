from MGPGtools.utils.common import run


def ogBuild(gfaFile, ogFile, threads):
    ODGIBuildCmd = (
        ["odgi", "build", "-g", gfaFile, "-s", "-o", ogFile]
        if threads == 1
        else ["odgi", "build", "-g", gfaFile, "-s", "-o", ogFile, "-t", threads]
    )
    run(ODGIBuildCmd)


def ogSort(ogFile, ogSortedFile, threads):
    ODGIBuildCmd = (
        ["odgi", "sort", "-i", ogFile, "-o", ogSortedFile]
        if threads == 1
        else ["odgi", "sort", "-i", ogFile, "-o", ogSortedFile, "-t", threads]
    )
    run(ODGIBuildCmd)


def ogExtract(ogFile, extractogFile, tPath, threads):
    ODGIExtractCmd = (
        [
            "odgi",
            "extract",
            "-i",
            ogFile,
            "-r",
            tPath,
            "-d",
            "3000",
            "-o",
            extractogFile,
        ]
        if threads == 1
        else [
            "odgi",
            "extract",
            "-i",
            ogFile,
            "-r",
            tPath,
            "-d",
            "3000",
            "-o",
            extractogFile,
            "-t",
            threads,
        ]
    )
    run(ODGIExtractCmd)


def ogExtractBed(ogFile, extractOgFile, bedFile, threads):
    ODGIExtractCmd = (
        [
            "odgi",
            "extract",
            "-i",
            ogFile,
            "-b",
            bedFile,
            "-d",
            "1000000",
            "-o",
            extractOgFile,
        ]
        if threads == 1
        else [
            "odgi",
            "extract",
            "-i",
            ogFile,
            "-b",
            bedFile,
            "-d",
            "1000000",
            "-o",
            extractOgFile,
            "-t",
            threads,
        ]
    )
    run(ODGIExtractCmd)


def ogPath(ogFile, threads):
    ODGIPathCmd = ["odgi", "paths", "-i", ogFile, "-H"] if threads == 1 else ["odgi", "paths", "-i", ogFile, "-H"]
    if_success, stdout, stderr = run(ODGIPathCmd)
    return if_success, stdout, stderr


def ogPathTsv(ogFile, tsvFile, threads):
    ODGIPathCmd = ["odgi", "paths", "-i", ogFile, "-H"] if threads == 1 else ["odgi", "paths", "-i", ogFile, "-H"]
    if_success, stdout, stderr = run(ODGIPathCmd)
    with open(tsvFile, "w") as f:
        f.write(stdout)


def ogView(ogFile, gfaFile, threads):
    ODGIViewCmd = (
        ["odgi", "view", "-i", ogFile, "-g"] if threads == 1 else ["odgi", "view", "-i", ogFile, "-g", "-t", threads]
    )
    if_success, stdout, stderr = run(ODGIViewCmd)
    with open(gfaFile, "w") as f:
        f.write(stdout)


def ogPosition(ogFile, positionFile, targetPath, threads):
    ODGIPositionCmd = (
        [
            "odgi",
            "position",
            "-i",
            ogFile,
            "-G",
            positionFile,
            "-r",
            targetPath,
        ]
        if threads == 1
        else [
            "odgi",
            "position",
            "-i",
            ogFile,
            "-G",
            positionFile,
            "-r",
            targetPath,
            "-t",
            threads,
        ]
    )
    if_success, stdout, stderr = run(ODGIPositionCmd)
    return if_success, stdout, stderr
