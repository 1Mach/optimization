import matplotlib.pyplot as plt
import numpy as np
import geatpy as ea
Cp_water=4200
V_inlet=89.4
T_inlet=251.35
Cp_air=1004.5
LWC=0.55
I_water=2500000
arealist=[]
betalist=[]
h_airlist=[]
pressurelist=[]
Twall_steadylist=[]
def read_original_data(filename):
	slist = []
	arealist = []
	datasetlist = []
	dataset2list = []
	dataset3list = []
	with open(filename, 'r') as f:
		lines = f.readlines()
		for data in lines:
			lines1 = data.strip('\n')
			datasetlist.append(lines1)
		for i in range(4, len(datasetlist)):
			dataset2list.append(datasetlist[i])
		for j in range(0, len(dataset2list)):
			dataset3list.append(dataset2list[j].split('\t'))
		for i in range(0, len(dataset3list) - 1):
			slist.append(float(dataset3list[i][0]))
			arealist.append(float(dataset3list[i][1]))
	return arealist

arealist.append(read_original_data('./original_data/area'))
area=np.asarray(arealist)
betalist.append(read_original_data('./original_data/beta'))
beta=np.asarray(betalist)
h_airlist.append(read_original_data('./original_data/heat_transfer_coefficient'))
h_air=np.asarray(h_airlist)
pressurelist.append(read_original_data('./original_data/pressure'))
pressure=np.asarray(pressurelist)
Twall_steadylist.append(read_original_data('./original_data/temperature'))
# print(Twall_steadylist)
Twall_steady=np.asarray(Twall_steadylist)
# print(Twall_steady[0,3])
Twall=np.ones((1,309))*309
# print(Twall)

def Q_extra(Twall, hair):
	T_inlet=251.35
	Qextra=hair*(Twall-T_inlet)
	return Qextra

def Q_out(Twall, Mout, area):
	Qout=Mout*Twall*Cp_water/area
	return Qout

def Q_in(Twall, Min, area):
	Qin=Min*Twall*Cp_water/area
	return Qin

def Q_imp(beta,area):
	Qimp=(V_inlet * beta * 0.55 * 0.001*area)*(V_inlet**2)/2+(V_inlet * beta * LWC * 0.001*area)*Cp_water*T_inlet
	Qimp=Qimp/area
	return Qimp

def M_imp(beta, area):
	Mimp=(V_inlet * beta * 0.55 * 0.001*area)
	return Mimp


def Q_aero(h_air):
	r=0.893131
	Qaero=r*(h_air*V_inlet**2)/(2*Cp_air)
	return Qaero

def M_evap(Twall, h_air, pressure, area):
	temp1 = 611 * pow(10, 7.45 * (Twall - 273.15) / (235 + Twall - 273.15))
	temp2=611 * pow(10, 7.45*(T_inlet- 273.15) / (235 + T_inlet - 273.15))
	hG = 0.622 *area* h_air / Cp_air
	Popera=0.0
	P_st=80004.023438
	Mevap=hG*(temp1 / ((pressure + Popera) - temp1) - temp2 / (P_st - temp2))
	return Mevap

def Q_evap(Mevap, Twall, area):
	Qevap=(Mevap*I_water+Mevap*Twall*Cp_water)/area
	return Qevap


def Q_cond(Qout, Qextra, Qevap, Qimp, Qaero, Qin):
	Qcond=Qout+Qextra+Qevap-Qimp-Qaero-Qin
	return Qcond

# def Water_film():
# 	Mevap=np.zeros((1, 309))
# 	Qevap=np.zeros((1, 309))
# 	for i in range(0, 309):
# 		Mevap[0,i]=M_evap(Twall[0,i], h_air[0,i], pressure[0,i], area[0,i])
# 		Qevap[0,i]=Q_evap(M_evap(Twall[0,i], h_air[0,i], pressure[0,i], area[0,i]), Twall[0,i], area[0,i])
# 	# print(Qevap)
# 	Mout=np.zeros((1, 309))
# 	Min=np.zeros((1, 309))
#
# 	Mout[0, 150]=M_imp(beta[0,150], area[0,150])-M_evap(Twall[0,150],h_air[0,150],pressure[0,150], area[0,150])
# 	for i in range(150, 308):
# 		Min[0,i+1]=Mout[0,i]
# 		Mout[0,i+1]=M_imp(beta[0,i+1], area[0,i])+Min[0,i+1]-M_evap(Twall[0,i+1],h_air[0,i+1],pressure[0,i+1], area[0,i+1])
# 		if Mout[0,i+1]<0:
# 			Mout[0,i+1]=0
#
# 	for i in range(0, 151):
# 		Min[0,150-i-1]=Mout[0,150-i]
# 		Mout[0,150-i-1]=M_imp(beta[0,150-i-1], area[0,i])+Min[0,150-i-1]-M_evap(Twall[0,150-i-1],h_air[0,150-i-1],pressure[0,150-i-1], area[0,150-i-1])
# 		if Mout[0,150-i-1]<0:
# 			Mout[0,150-i-1]=0
# 	return Mout, Min
#





def aim(Twall,Mout, Min, area, h_air, pressure, beta ):  # 传入种群染色体矩阵解码后的基因表现型矩阵
	return Q_out(Twall, Mout, area)+Q_extra(Twall, h_air)+Q_evap(M_evap(Twall, h_air, pressure, area),Twall, area)-Q_imp(beta, area)-Q_aero(h_air)-Q_in(Twall, Min, area)


class MyProblem(ea.Problem):
	def __init__(self):
		name='MyProblem'
		M=1
		maxormins=[1]
		Dim=309
		varTypes=[0]*Dim
		lb=[290 for x in range(0, 309)]
		ub=[320 for x in range(0, 309)]
		lbin=[0 for x in range(0,309)]
		ubin=[1 for x in range(0, 309)]
		ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
	def aimFunc(self, pop):#目标函数
		aim2=0
		Mevap = np.zeros((1, 309))
		Qevap = np.zeros((1, 309))
		Mout = np.zeros((1, 309))
		Min = np.zeros((1, 309))

		Vars=pop.Phen#种群表现型
		x = np.zeros((16, 309))#309cells数， 16种群规模
		Min=np.zeros((16,309))
		Mout=np.zeros((16,309))
		for i in range(0,309):
			x[:,[i]] = Vars[:, [i]]
		for j in range(0, 309):
			Mevap[0, j] = M_evap(x[-1,j], h_air[0, j], pressure[0, j], area[0, j])
			Qevap[0, j] = Q_evap(M_evap(x[-1,j], h_air[0, j], pressure[0, j], area[0, j]), x[-1,j], area[0, j])
		# print(x[-1,[307]])
		Mout[:, [150]] = M_imp(beta[0, 150], area[0, 150]) - M_evap(x[:, [150]], h_air[0, 150], pressure[0, 150], area[0, 150])
		# print(Mout[0, 150])
		# print(x[-1, 150])
		for i in range(150, 308):#mout和min没有进行进化，只是利用了种群中的一个个体，这是有问题滴
			Min[:, [i + 1]] = Mout[:, [i]]
			Mout[:, [i + 1]] = M_imp(beta[0, i + 1], area[0, i+1]) + Min[:,[ i + 1]] - M_evap(x[-1, [i+1]], h_air[0, i + 1],pressure[0, i + 1],	area[0, i + 1])
			if Mout[:, [i + 1]].all() < 0:
				Mout[:, [i + 1]] = 0

		for i in range(0, 151):
			Min[:, [150 - i - 1]] = Mout[:, [150 - i]]
			Mout[:, [150 - i - 1]] = M_imp(beta[0, 150 - i - 1], area[0, i]) + Min[:, [150 - i - 1]] - M_evap(x[:,[150-i-1]], h_air[0, 150 - i - 1], pressure[0, 150 - i - 1], area[0, 150 - i - 1])
			if Mout[:, [150 - i - 1]].all() < 0:
				Mout[:, [150 - i - 1]] = 0



		for i in range(0, 309):
			aim2+=aim(x[:, [i]], Mout[:,[i]], Min[:,[i]], area[0,i],h_air[0,i], pressure[0,i], beta[0,i])*area[0, i]
			# print(aim(x[:, [i]], Mout[0,i], Min[0,i], area[0,i],h_air[0,i], pressure[0,i], beta[0,i]))
		pop.ObjV=aim2
		pop.CV=np.hstack([
			              # M_imp(beta[0, 306 ], area[0, 306]) + Mout[0, 305] - M_evap(x[:, [306]], h_air[0, 306], pressure[0, 306],  area[0, 306]),
						  # M_imp(beta[0, 307 ], area[0, 307]) + Mout[0, 306] - M_evap(x[:, [307]], h_air[0, 307],pressure[0, 307],	area[0, 307]),
						  M_imp(beta[0, 308 ], area[0, 308]) + Mout[:, [307]] - M_evap(x[:, [308]], h_air[0, 308],pressure[0, 308],	area[0, 308]),
						  M_imp(beta[0, 0], area[0, 0]) + Mout[:, [1]] - M_evap(x[:, [0]], h_air[0, 0], pressure[0, 0],   area[0, 0]),
						  # M_imp(beta[0, 1], area[0, 1]) + Mout[0, 2] - M_evap(x[:, [1]], h_air[0, 1], pressure[0, 1],   area[0, 1]),
						  # M_imp(beta[0, 2], area[0, 2]) + Mout[0, 3] - M_evap(x[:, [2]], h_air[0, 2], pressure[0, 2],   area[0, 2]),
						  # np.abs(x[:, [306]] - 300),
						  np.abs(x[:,[307]]-x[:,[308]])-5
						  # np.abs(x[:,[307]]-300),
						  # np.abs(x[:,[308]]-380)
						  ])
		# print( M_evap(x[:, [308]], h_air[0, 308],pressure[0, 308],area[0, 308]))
		# print( M_imp(beta[0, 308 ], area[0, 308]))
		# print(x[-1, 308])
		# print(M_imp(beta[0, 150], area[0, 150]))
		# print(M_evap(x[-1, 150], h_air[0, 150], pressure[0, 150],  area[0, 150]))
		print(Mout[-1,:])
		# y=np.arange(0, 309).reshape(1, 309)
		# plt.ion()
		# plt.plot(y, Mout)




		# pop.CV = np.hstack([np.abs(M_evap(x[:, [308]], h_air[0, 308],pressure[0, 308],	area[0, 308]))])

