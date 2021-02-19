import numpy as np

for i in range(0, 20):
	print('-(M_imp(beta[:, [%d]], area[:, [%d]]) + Mout[:, [%d]] - M_evap(x[:, [%d]], h_air[:, [%d]], pressure[:, [%d]],area[:, [%d]])),' %(308-i, 308-i,307-i, 308-i,308-i, 308-i,308-i))
	print('M_imp(beta[:, [%d] ], area[:, [%d]]) + Mout[:, [%d]] - M_evap(x[:, [%d]], h_air[:, [%d]],pressure[:, [%d]],	area[:, [%d]])-0.00002,' %(308-i, 308-i,307-i, 308-i,308-i, 308-i,308-i))