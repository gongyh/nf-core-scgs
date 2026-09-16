########################### f1,nmi等指标计算 ##########################
#################################################################
from sklearn import metrics
from sklearn.metrics import silhouette_score,davies_bouldin_score,normalized_mutual_info_score,adjusted_rand_score,fowlkes_mallows_score
#打开txt文件并读取其中的内容
with open("/media/ubuntu/abc/csm/ablabtionLabs/cami_low/final_files/camilow_final_prediction.txt", 'r') as f: # 预测标签
    content = f.read()

with open("/media/ubuntu/abc/csm/ablabtionLabs/cami_low/camilowover2k_40gs_true.txt", 'r') as fd: # 真实标签
    contentd = fd.read()
# 将内容按行拆分，并保存为列表
label_pred = content.splitlines() # 预测标签
label_true = contentd.splitlines() # 真实标签
accuracy=metrics.accuracy_score(label_true, label_pred)
f1=metrics.f1_score(
          label_true, label_pred, average="macro", zero_division=0
      )
precision=metrics.precision_score(
            label_true, label_pred, average="macro", zero_division=0
        )
recall=metrics.recall_score(
            label_true, label_pred, average="macro", zero_division=0
        )

ari = adjusted_rand_score(label_true, label_pred)
nmi = normalized_mutual_info_score(label_true, label_pred)
fmi = fowlkes_mallows_score(label_true, label_pred)
print("*** Pred results ***\n""acc:%f\nf1:%f\npre:%f\nrecall:%f\nari:%f\nnmi:%f\nfmi:%f\n" %(accuracy,f1,precision,recall,ari,nmi,fmi))