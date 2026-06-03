import parsecontigs
import vambtools
import parsebam
import argparse
from typing import Optional, Tuple, Union, cast
from pathlib import Path
import numpy as np
import sys
from loguru import logger
from math import isfinite
import time
import os
import pycoverm

parser = argparse.ArgumentParser(description="提取tnf,rpkm")
parser.add_argument('-od','--output_dir',type=str,help="输出路径")
parser.add_argument('-fd','--fasta_file',type=str,help="fasta输入路径")
parser.add_argument('-bam','--bam_file',type=str,help="bam输入路径")
args = parser.parse_args()
base_name = os.path.splitext(os.path.basename(args.fasta_file))[0]
# 计算fasta文件的tnf向量（4mer频率降维103）输出为(nx103的ndarry),n是contigs的个数
class FASTAPath(type(Path())):
    __slots__ = []
    pass
class CompositionPath(type(Path())):
    __slots__ = []
    pass
class CompositionOptions:
    __slots__ = ["path", "min_contig_length"]

    def __init__(
        self,
        fastapath: Optional[Path],
        npzpath: Optional[Path],
        min_contig_length: int,
    ):
        assert isinstance(fastapath, (Path, type(None)))
        assert isinstance(npzpath, (Path, type(None)))
        assert isinstance(min_contig_length, int)

        if min_contig_length < 250:
            raise argparse.ArgumentTypeError(
                "Minimum contig length must be at least 250"
            )

        if not (fastapath is None) ^ (npzpath is None):
            raise argparse.ArgumentTypeError(
                "Must specify either FASTA or composition path"
            )

        for path in (fastapath, npzpath):
            if path is not None and not path.is_file():
                raise FileNotFoundError(path)

        if fastapath is not None:
            self.path = FASTAPath(fastapath)
        else:
            assert npzpath is not None
            self.path = CompositionPath(npzpath)
        self.min_contig_length = min_contig_length

# fastapath = Path("/media/ubuntu/abc/csm/KELIN/My_DnabertTwo/myvae/testdata/sample1.fna")
# 如果事先传入kernal.npz文件（一种用于生成tnf的约束），则使用文件，如果没有则重新写入composition.npz
# npzpath = Path("/media/ubuntu/abc/csm/KELIN/My_DnabertTwo/myvae/mainfiles/kernel.npz")
# options = CompositionOptions(fastapath,npzpath,2000) 
# options = CompositionOptions(fastapath,None,2000) #  传入两个参数，一个是fasta文件路径，一个是min_contig_length(最短长度)使用2000代替options.min_contig_length

with vambtools.Reader(args.fasta_file) as file:
    # composition = parsecontigs.Composition.from_file(
    #             file, minlength=options.min_contig_length
    #         )
    composition = parsecontigs.Composition.from_file(
                file, minlength=4
            )


composition.save(Path(args.output_dir).joinpath(f"{base_name}_tnf.npz"))

binsplitter = vambtools.BinSplitter.inert_splitter()
binsplitter.initialize(composition.metadata.identifiers)

logger.remove()  # 移除默认的 handler
logger.add(sys.stdout, format="{level} | {message}")  # 只保留等级和消息
# 如果观察到有任何 contigs 的长度小于设定的阈值时，作出提醒
if not np.all(composition.metadata.mask):
        n_removed = len(composition.metadata.mask) - np.sum(composition.metadata.mask)
        message = (
            f"The minimum sequence length has been set to {2000}, "
            f"but {n_removed} sequences fell below this threshold and was filtered away."
            "\nBetter results are obtained if the sequence file is filtered to the minimum "
            "sequence length before mapping.\n"
        )
        logger.opt(raw=True).info("\n")
        logger.warning(message)
begintime = time.time()
elapsed = round(time.time() - begintime, 2)
logger.info(
    f"\tKept {composition.count_bases()} bases in {composition.nseqs} sequences"
)
logger.info(f"\tProcessed TNF in {elapsed} seconds.\n")
print("TNF信息如下：")
print(composition.matrix)
print(composition.matrix.shape) # (2510 * 103)

########计算rpkm需要利用TNF的信息#####################################################################################
class AbundancePath(type(Path())):
    pass
class AbundanceOptions:
    __slots__ = ["path", "min_alignment_id", "refcheck"]

    def __init__(
        self,
        bampaths: Optional[list[Path]],
        abundancepath: Optional[Path],
        min_alignment_id: Optional[float],
        refcheck: bool,
    ):
        assert isinstance(bampaths, (list, type(None)))
        assert isinstance(abundancepath, (Path, type(None)))
        assert isinstance(min_alignment_id, (float, type(None)))
        assert isinstance(refcheck, bool)

        # Make sure only one RPKM input is there
        if not (bampaths is not None) + (abundancepath is not None) == 1:
            raise argparse.ArgumentTypeError(
                "Must specify exactly one of BAM files or abundance NPZ file input"
            )

        if abundancepath is not None:
            if not abundancepath.is_file():
                raise FileNotFoundError(
                    f'Not an existing non-directory file: "{str(abundancepath)}"'
                )
            self.path = AbundancePath(abundancepath)

        elif bampaths is not None:
            for bampath in bampaths:
                if not bampath.is_file():
                    raise FileNotFoundError(
                        f'Not an existing non-directory file: "{str(bampath)}"'
                    )
                if not pycoverm.is_bam_sorted(str(bampath)):
                    raise ValueError(f"Path {bampath} is not sorted by reference.")
            self.path = bampaths

        if min_alignment_id is not None:
            if bampaths is None:
                raise argparse.ArgumentTypeError(
                    "If minid is set, RPKM must be passed as bam files"
                )
            if (
                not isfinite(min_alignment_id)
                or min_alignment_id < 0.0
                or min_alignment_id > 1.0
            ):
                raise argparse.ArgumentTypeError(
                    "Minimum nucleotide ID must be in [0,1]"
                )
            self.min_alignment_id = min_alignment_id
        else:
            self.min_alignment_id = 0.0

        self.refcheck = refcheck
begintime = time.time()

bampaths = [Path(args.bam_file)] # Optional[list[Path]] # bam文件路径
min_alignment_id = 0.0
# refcheck = True # 如果refcheck为true,则打印参考哈希值，否则打印none

abundance_options = AbundanceOptions(bampaths,None,0.0,True)

path = abundance_options.path
comp_metadata = composition.metadata # 数据类型为parsecontigs.CompositionMetaData，由tnf:compositon.metadata给出。
nthreads = 64
abundance = parsebam.Abundance.from_files(
            path,
            Path(args.output_dir).joinpath("tmp").joinpath("pycoverm"),
            comp_metadata,
            abundance_options.refcheck,
            abundance_options.min_alignment_id,
            nthreads,
        )
abundance.save(Path(args.output_dir).joinpath("rpkm.npz"))

        
for i, samplename in enumerate(abundance.samplenames):
    logger.info(f"\t{i:>6}: {samplename}")

