#!/bin/bash

#SBATCH --job-name=fur
#SBATCH --output=fur_output.txt
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --partition=gpujl
#SBATCH --gres=gpu:4
# =====================
# 运行环境说明：
# 脚本功能：运行从spades组装、特征提取、VAE融合、聚类、checkm2评估的全流程
# 参数顺序：<input_dir> <output_dir>
# =====================
# if ! command -v conda &>/dev/null; then
#     echo "ERROR: Conda not found!" >&2
#     exit 1
# fi
# eval "$(conda shell.bash hook 2>/dev/null)" || {
#     conda init bash
#     source ~/.bashrc
#     eval "$(conda shell.bash hook)"
# }

# ===== 环境清理部分 =====
# 清除所有可能冲突的 conda 环境变量
unset CONDA_PREFIX
unset CONDA_DEFAULT_ENV
unset CONDA_PROMPT_MODIFIER
unset CONDA_EXE
unset CONDA_PYTHON_EXE
unset CONDA_SHLVL

# 清除 Perl 相关环境变量
unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
unset PERL_MB_OPT
unset PERL_MM_OPT

# 清理 PATH 中的 conda 路径
export PATH=$(echo "$PATH" | tr ':' '\n' | grep -v "/home/wangjingyuan/anaconda3" | tr '\n' ':' | sed 's/:$//')
# ===== 环境清理结束 =====

clean_env() {
    unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_PROMPT_MODIFIER
    unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
    export PATH=$(echo "$PATH" | tr ':' '\n' | grep -v "/home/wangjingyuan/anaconda3" | tr '\n' ':' | sed 's/:$//')
}


# 使用正确的 conda 路径重新初始化
#CONDA_PATH="/home/wangjingyuan/anaconda3"

# 提示用户输入 conda 路径
echo "请输入 conda 初始化路径："
read CONDA_PATH

# 使用用户输入的路径
echo "使用 conda 路径: $CONDA_PATH"

if [ -f "$CONDA_PATH/etc/profile.d/conda.sh" ]; then
    source "$CONDA_PATH/etc/profile.d/conda.sh"
else
    export PATH="$CONDA_PATH/bin:$PATH"
fi





start_time=$(date +%s)

# 获取脚本所在目录并切换到项目根目录
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
cd "$PROJECT_ROOT" || exit 1

# 检查参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir>  <output_dir>"
    exit 1
fi

# 输入参数
input_dir=$1
output_dir=$2

echo "$PROJECT_ROOT"

conda activate vambnew
# Step 1: 运行组装脚本

python $PROJECT_ROOT/scripts/spades.py \
   -i "$input_dir" -t 66

# Step 2: 自动查找fasta和bam
fasta_file=$(find "$input_dir" -type f -name "*over2k.fasta" | head -n 1)
bam_file=$(find "$input_dir" -type f -name "*_sorted.bam" | head -n 1)

if [ ! -f "$fasta_file" ] || [ ! -f "$bam_file" ]; then
    echo "ERROR: Required .fasta or .bam file not found in $input_dir"
    exit 1
fi

# 提取文件基本名和父目录
BASE_NAME=$(basename "$fasta_file" .fasta)
PARENT_DIR=$(dirname "$fasta_file")

#生成预训练数据集
# 创建输出文件夹
mkdir -p "$output_dir"
# 创建onlyseq序列文件输出路径
seq_file="${output_dir}/${BASE_NAME}_over2kseq.txt"
# 创建fpf特征输出路径
fpf_file="${output_dir}/${BASE_NAME}_fpf.npy"

conda deactivate

conda activate dnaberts

pretrain_dir="${output_dir}/furpretrainfiles"
mkdir -p "$pretrain_dir"

python $PROJECT_ROOT/scripts/generate_pretrain_data.py \
    -fasta "$fasta_file" -outdir "$pretrain_dir"
echo "Pretraining data saved in: ${pretrain_dir}"

export PATH_TO_DATA_DICT=$pretrain_dir
export TRAIN_FILE=train_set.csv
export VAL_FILE=val_set.csv
export model_dir=$pretrain_dir/model_output

echo "$fasta_file"
echo "$pretrain_dir"



#启动训练脚本
cd "$PROJECT_ROOT/DNABERT_S-main/train/pretrain"
export NUMEXPR_MAX_THREADS=32
python main.py \
  --resdir ${model_dir} \
  --datapath ${PATH_TO_DATA_DICT} \
  --train_dataname ${TRAIN_FILE} \
  --val_dataname ${VAL_FILE} \
  --seed 1 \
  --logging_step 2000 \
  --logging_num 15 \
  --max_length 2000 \
  --train_batch_size 16 \
  --val_batch_size 128 \
  --lr 1e-06 \
  --lr_scale 100 \
  --epochs 3 \
  --feat_dim 128 \
  --temperature 0.05 \
  --con_method same_species \
  --mix \
  --mix_alpha 1.0 \
  --mix_layer_num -1 \
  --curriculum

# 激活dnaberts环境进行特征提取
#__conda_setup="$('/home/ubuntu/mambaforge/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
export CONDA_ROOT="${CONDA_ROOT:-$HOME/mambaforge}"  # 默认值: $HOME/mambaforge
__conda_setup="$("${CONDA_ROOT}/bin/conda" 'shell.bash' 'hook' 2> /dev/null)"

if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/ubuntu/mambaforge/etc/profile.d/conda.sh" ]; then
        . "/home/ubuntu/mambaforge/etc/profile.d/conda.sh"
    else
        export PATH="/home/ubuntu/mambaforge/bin:$PATH"
    fi
fi
unset __conda_setup
#conda activate dnaberts

# 1.运行特征提取文件
python $PROJECT_ROOT/scripts/featureExtract_gpu_2.py -md "$model_dir" -fd "$fasta_file" -sd "$seq_file" -dd "$fpf_file" 
echo "fpf_features output: ${fpf_file}"
conda deactivate

# 2.激活vamb环境进行fpf，rpkm提取
conda activate vambnew
python $PROJECT_ROOT/myvae/mainfiles/calc_tnf_and_rpkm_2.py -od "$output_dir" -fd "$fasta_file" -bam "$bam_file"

# 将会在output_dir文件夹下生成tnf.npz和abundance.npz文件
tnf_file="${output_dir}/${BASE_NAME}_tnf.npz"
rpkm_file="${output_dir}/rpkm.npz"
echo "tnf_features output: ${tnf_file}"

# 3.特征融合
# 创建vae特征输出路径
vaef_file="${output_dir}/${BASE_NAME}_vae_features.npy"
python $PROJECT_ROOT/myvae/mainfiles/vaeTest_2.py -dd "$fpf_file" -td "$tnf_file" -rd "$rpkm_file" -vd "$vaef_file"
echo "vaefeatures output: ${vaef_file}"

# 4.计算初始聚类个数
kmer_file="${output_dir}/4mer.csv"
seqid_file="${output_dir}/seqid.csv"
python $PROJECT_ROOT/scripts/calculate_kmer_multi_thread_2.py "$fasta_file" "$kmer_file" "$seqid_file" 4
conda deactivate

conda activate copygen
cvf_file="${output_dir}/cluster_value"

python $PROJECT_ROOT/marker_gene/src/marker_gene_utils.py -kf "$kmer_file" -cf "$fasta_file" -cvf "$cvf_file"
echo "cluster_value output：${cvf_file}"
conda deactivate

# 5.聚类
# 创建预测标签输出路径
LABEL_FILE="${output_dir}/${BASE_NAME}_prinum.txt"
# 事先创建分箱结果文件夹
BINS_OUTPUT_DIR="${output_dir}/${BASE_NAME}_bins"
mkdir -p "$BINS_OUTPUT_DIR"
conda activate dnaberts
python $PROJECT_ROOT/scripts/myCluster_2.py -vd "$vaef_file" -ld "$LABEL_FILE" -fd "$fasta_file" -bd "$BINS_OUTPUT_DIR" -cvf "$cvf_file"
echo "bins output: ${BINS_OUTPUT_DIR}"
conda deactivate

# 6.运行checkm2指标
conda activate checkm2
checkm2 predict -t 64 -x fasta --input "${BINS_OUTPUT_DIR}" --output-directory "${output_dir}/${BASE_NAME}_checkm2"
conda deactivate

# 7.统计箱子质量
# 定义 .tsv 文件路径
tsv_filea="${output_dir}/${BASE_NAME}_checkm2/quality_report.tsv"
# 统计 Completeness >= 90 并且 Contamination < 5 的数量
count1a=$(awk -F '\t' 'NR > 1 && $2 >= 90 && $3 < 5 {count++} END {print count+0}' "$tsv_filea")
# 统计 Completeness >= 80 且 Contamination < 5 的数量
count2a=$(awk -F '\t' 'NR > 1 && $2 >= 80 && $3 < 5 {count++} END {print count+0}' "$tsv_filea")
# 统计 Completeness >= 70 且 Contamination < 5 的数量
count3a=$(awk -F '\t' 'NR > 1 && $2 >= 70 && $3 < 5 {count++} END {print count+0}' "$tsv_filea")
# 统计 Completeness >= 60 且 Contamination < 5 的数量
count4a=$(awk -F '\t' 'NR > 1 && $2 >= 60 && $3 < 5 {count++} END {print count+0}' "$tsv_filea")
# 统计 Completeness >= 50 且 Contamination < 5 的数量
count5a=$(awk -F '\t' 'NR > 1 && $2 >= 50 && $3 < 5 {count++} END {print count+0}' "$tsv_filea")

# 统计mags标准下的
# 统计 Completeness > 90 且 Contamination < 5 的数量
count6a=$(awk -F '\t' 'NR > 1 && $2 > 90 && $3 < 5 {count++} END {print count+0}' "$tsv_filea")
# 统计 Completeness >= 50 且 Contamination < 10 的数量
count7a=$(awk -F '\t' 'NR > 1 && $2 >= 50 && $3 < 10 {count++} END {print count+0}' "$tsv_filea")

echo "Completeness >= 90 and Contamination < 5: $count1a"
echo "Completeness >= 80 and Contamination < 5: $count2a"
echo "Completeness >= 70 and Contamination < 5: $count3a"
echo "Completeness >= 60 and Contamination < 5: $count4a"
echo "Completeness >= 50 and Contamination < 5: $count5a"

echo "MAGs标准下的指标结果："
echo "Completeness > 90 and Contamination < 5: $count6a"
echo "Completeness >= 50 and Contamination < 10: $count7a"

# 记录脚本结束时间
end_time=$(date +%s)
# 计算总运行时间
runtime=$((end_time - start_time))
echo "运行时间: $runtime seconds"