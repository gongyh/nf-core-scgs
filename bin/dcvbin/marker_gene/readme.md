
运行前切换到"copygen"环境下(test环境配置：step1:conda create -n test python=3.6.13 step2:conda install scikit-learn==0.22.2.post1 )
- 运行方法：

  ```python
  yanziming@server1:~/vicent/contrast_experiment/marker_gene$ python ./src/marker_gene_utils.py
  ```

- 需要两个文件（在marker_gene_utils.py文件中设置）：

  fasta文件：'/home/yanziming/csm/data_set/sharon/unknownseq.fasta'

  该fasta文件的kmer频率文件：'/home/yanziming/csm/data_set/sharon/kmer.csv'


- 运行会报错：UnboundLocalError: local variable 'candK' referenced before assignment
  需要修改marker_gene中一些文件的执行权限，可在.err文件（报错文件）中查看。（文件路径和传入的fasta文件相同）
 1.auxiliary/hmmer-3.3/src/hmmsearch
 2.auxiliary/FragGeneScan1.31/run_FragGeneScan.pl
 3.auxiliary/FragGeneScan1.31/FragGeneScan
 4.auxiliary/test_getmarker.pl


- 4-mer频率归一化也会影响最后的结果

  ```python
  x_contigs = normalize(x_contigs, norm='l2', axis=1) # 会影响最终的结果，全部序列归一化后结果为12，归一化之前是14。未知序列归一化后结果为4，归一化之前也为4。
  ```

  全部序列单拷贝结果归一化后结果为**12**个，归一化之前是**14**个，未知序列归一化后结果为**4**，归一化之前也为**4**

