import numpy as np
import matplotlib.pyplot as plt
import re
#Read related file
filename='./original_data/beta'
s=[]
beta=[]
dataset=[]
dataset2=[]
dataset3=[]
with open(filename, 'r') as f:
	lines=f.readlines()
	for data in lines:
		lines1=data.strip('\n')
		lines2=data.split('\t')
		dataset.append(lines1)
	for i in range(4, len(dataset)):
		dataset2.append(dataset[i])
	for j in range(0, len(dataset2)):
		dataset3.append(dataset2[j].split('\t'))
	for i in range(0, len(dataset3)-1):
		s.append(float(dataset3[i][0]))
		beta.append(float(dataset3[i][1]))
	print(s)
	print(beta)

	
plt.plot(s, beta,'ro')
plt.show()

