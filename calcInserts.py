import numpy as np 

# col1: occurence; col2: read length
L001 = np.loadtxt('../data/GEO_Submission/AL1_S1_L001_R2_001.length.txt', dtype='float')
print(np.mean(L001[:,1]), np.std(L001[:,1]))

L002 = np.loadtxt('../data/GEO_Submission/AL1_S1_L002_R2_001.length.txt', dtype='float')
print(np.mean(L002[:,1]), np.std(L002[:,1]))