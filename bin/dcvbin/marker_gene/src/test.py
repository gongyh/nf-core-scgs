import numpy as np

data = np.loadtxt('/home/yanziming/vicent/data_set/healthy_human_fecal/shotgun/SRR5298275/spades_out/4mer.csv',
                  delimiter=',', usecols=range(1, 137))
print(data)
