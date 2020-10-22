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

def Water_film():
	Mevap=np.zeros((1, 309))
	Qevap=np.zeros((1, 309))
	for i in range(0, 309):
		Mevap[0,i]=M_evap(Twall[0,i], h_air[0,i], pressure[0,i], area[0,i])
		Qevap[0,i]=Q_evap(M_evap(Twall[0,i], h_air[0,i], pressure[0,i], area[0,i]), Twall[0,i], area[0,i])
	# print(Qevap)
	Mout=np.zeros((1, 309))
	Min=np.zeros((1, 309))

	Mout[0, 150]=M_imp(beta[0,150], area[0,150])-M_evap(Twall[0,150],h_air[0,150],pressure[0,150], area[0,150])
	for i in range(150, 308):
		Min[0,i+1]=Mout[0,i]
		Mout[0,i+1]=M_imp(beta[0,i+1], area[0,i])+Min[0,i+1]-M_evap(Twall[0,i+1],h_air[0,i+1],pressure[0,i+1], area[0,i+1])
		if Mout[0,i+1]<0:
			Mout[0,i+1]=0

	for i in range(0, 151):
		Min[0,150-i-1]=Mout[0,150-i]
		Mout[0,150-i-1]=M_imp(beta[0,150-i-1], area[0,i])+Min[0,150-i-1]-M_evap(Twall[0,150-i-1],h_air[0,150-i-1],pressure[0,150-i-1], area[0,150-i-1])
		if Mout[0,150-i-1]<0:
			Mout[0,150-i-1]=0
	return Mout, Min

Mout, Min=Water_film()
# print(Mout)
# y=np.arange(0, 309, 1)
# plt.plot(y, Mout[0])
# plt.show()


def aim(Twall):  # 传入种群染色体矩阵解码后的基因表现型矩阵
	return Q_out(Twall, Mout[0,10], area[0,10])+Q_extra(Twall, h_air[0,10])+Q_evap(M_evap(Twall, h_air[0,10], pressure[0,10], area[0,10]),Twall, area[0,10])-Q_imp(beta[0,10], area[0,10])-Q_aero(h_air[0,10])-Q_in(Twall, Min[0,10], area[0,10])



class MyProblem(ea.Problem):
	def __init__(self):
		name='MyProblem'
		M=1
		maxormins=[1]
		Dim=309
		varTypes=[0]*Dim
		lb=[290 for x in range(0, 309)]
		ub=[309 for x in range(0, 309)]
		lbin=[0 for x in range(0,309)]
		ubin=[1 for x in range(0, 309)]
		ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
	def aimFunc(self, pop):#目标函数
		aim2=0
		Vars=pop.Phen
		x={}
		for i in range(0,309):
			x[i] = Vars[:, [i]]
			aim2+=aim(x[i])
		pop.ObjV=aim2
