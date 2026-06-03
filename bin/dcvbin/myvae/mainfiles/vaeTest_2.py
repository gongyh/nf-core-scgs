from myencode import VAE, make_dataloader
import numpy as np
from pathlib import Path
import argparse
parser = argparse.ArgumentParser(description="用于提取VAE特征:fpf_dir、tnf_dir、rpkm_dir、vaefeatures_dir")
parser.add_argument('-dd','--fpf_file',type=str,help="指定fpf特征路径")
parser.add_argument('-td','--tnf_file',type=str,help="指定tnf特征路径")
parser.add_argument('-rd','--rpkm_file',type=str,help="指定rpkm特征路径")
parser.add_argument('-vd','--vaef_file',type=str,help="指定vae特征输出路径")
args = parser.parse_args()
# 实例化 VAE 对象
vae = VAE(nsamples=1)

# 创建数据加载器
# 假设 depths、tnf 和 lengths 是你的数据，你需要定义 make_dataloader 函数来生成一个数据加载器
# tnf是4-mer频率256维通过降维方式映射成103维
############################################################################################
# 构建临时的dnabert模型向量DBF(dna bert features)
# fpf = np.random.random((2510, 768)).astype(np.float32) # 2510 x 768
fpf = np.load(args.fpf_file,allow_pickle=True) # 2510 x 768

############################################################################################

# 加载tnf
tnf = np.load(args.tnf_file,allow_pickle=True) # 2510 x 103
# array_names = tnf.files # ['matrix', 'identifiers', 'lengths', 'mask', 'minlength']

# 加载rpkm
rpkm = np.load(args.rpkm_file,allow_pickle=True) # 2510 x 1
# rpkm_array_names = rpkm.files # ['matrix', 'samplenames', 'minid', 'refhash']

# 长度length
lens =  tnf['lengths'] # [2012 2212 3515 .... 25255] 

# rpkm: RPKM matrix (N_contigs x N_samples)
# tnf: TNF matrix (N_contigs x N_TNF)
# dbf: TBF matrix (N_contigs x 768)
# lengths: Numpy array of sequence length (N_contigs)
dataloader = make_dataloader(rpkm['matrix'], tnf['matrix'], fpf, lens)

# 训练 VAE 模型
# modelsave_path = Path("/media/ubuntu/abc/csm/KELIN/MY_DNABERTS/mydata/SRR952674/vaeoutput/modelover2k_dbs.pt")
vae.trainmodel(dataloader)

# 对数据进行编码，得到潜在表示
latent = vae.encode(dataloader)
np.save(args.vaef_file,latent)


# 打印潜在表示的形状
print(latent.shape)  # 输出：(2510, 32)

