import numpy as np
import math as m
NIND=100
n=20
T0=np.ones((NIND, n))*273.15
T_s=np.zeros((NIND, n))
T_wall=np.zeros((NIND, n))
i=10
# for j in range(0, NIND):
#     print(m.exp(T0[j, i]))
T_array=[]
T_array=T0[0]
#print(abs(T_s[:, [i]]-T0[:, [0]]))
print((T_array<=274).all())