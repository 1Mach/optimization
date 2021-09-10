

import matplotlib.pyplot as plt
import numpy as np
import geatpy as ea
Cp_water=4200
V_inlet=89.4
T_inlet=251.35
Cp_air=1004.5
LWC=0.55
NIND=1
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
area=np.repeat(area, NIND, axis=0)#复制行
# print(area)
betalist.append(read_original_data('./original_data/beta'))
beta=np.asarray(betalist)
beta=np.repeat(beta, NIND, axis=0)
h_airlist.append(read_original_data('./original_data/heat_transfer_coefficient'))
h_air=np.asarray(h_airlist)
h_air=np.repeat(h_air, NIND, axis=0)
pressurelist.append(read_original_data('./original_data/pressure'))
pressure=np.asarray(pressurelist)
pressure=np.repeat(pressure, NIND, axis=0)


Twall_steadylist.append(read_original_data('./original_data/temperature'))
# print(Twall_steadylist)
Twall_steady=np.asarray(Twall_steadylist)#计算的稳态的温度
# print(Twall_steady[0,3])
Twall=np.ones((1,309))*306#转换成矩阵，初始温度
# print(Twall)

def Q_cond(q1, q2, q3, q4, q5, q6, q7):#First, the partition is random(the cell number in every part is random)
	Qcond = np.zeros((1, 309))
	for i in range(0, 40):
		Qcond[0, i] = q1
	for i in range(40,82):
		Qcond[0, i]=q2
	for i in range(82, 135):
		Qcond[0, i]=q3
	for i in range(135, 186):
		Qcond[0, i]=q4
	for i in range(186, 230):
		Qcond[0, i]=q5
	for i in range(230, 274):
		Qcond[0, i]=q6
	for i in range(274, 309):
		Qcond[0, i]=q7
	return Qcond


def Q_extra(Twall, hair):#计算对流换热损失
	T_inlet=251.35
	Qextra=hair*(Twall-T_inlet)
	return Qextra

def Q_out(Twall, Mout, area):#计算流出热流密度
	Qout=Mout*Twall*Cp_water
	return Qout

def Q_in(Twall, Min, area):#计算流入热流密度
	Qin=Min*Twall*Cp_water
	return Qin

def Q_imp(beta,area):#计算撞击热流密度
	Qimp=(V_inlet * beta * 0.55 * 0.001*area)*(V_inlet**2)/2+(V_inlet * beta * LWC * 0.001*area)*Cp_water*T_inlet
	return Qimp

def M_imp(beta, area):#计算撞击质量
	Mimp=V_inlet * beta * LWC * 0.001*area
	return Mimp


def Q_aero(h_air):#计算气动加热热流密度
	r=0.893131
	Qaero=r*(h_air*V_inlet**2)/(2*Cp_air)
	return Qaero

def M_evap(Twall, h_air, pressure, area):#计算蒸发质量
	temp1 = 611 * pow(10, 7.45 * (Twall - 273.15) / (235 + Twall - 273.15))
	temp2=611 * pow(10, 7.45*(T_inlet- 273.15) / (235 + T_inlet - 273.15))
	hG = 0.622 *area* h_air / Cp_air
	Popera=0.0
	P_st=80004.023438
	Mevap=hG*(temp1 / ((pressure + Popera) - temp1) - temp2 / (P_st - temp2))
	return Mevap

def Q_evap(Mevap, Twall, area):#计算蒸发热流密度
	Qevap=Mevap*I_water+Mevap*Twall*Cp_water
	return Qevap



def Water_film():  #计算水膜
	Mevap=np.zeros((0, 309))   #初始化
	Qevap=np.zeros((0, 309))    #初始化
	for i in range(0, 309):
		Mevap[:,[i]]=M_evap(Twall[:,[i]], h_air[:,[i]], pressure[:,[i]], area[:,[i]])   #计算蒸发质量
		Qevap[:,[i]]=Q_evap(M_evap(Twall[:,[i]], h_air[:,[i]], pressure[:,[i]], area[:,[i]]), Twall[:,[i]], area[:,[i]])   #计算蒸发能量
	# print(Qevap)
	Mout=np.zeros((NIND, 309))  #初始化
	Min=np.zeros((NIND, 309))    #初始化

	Mout[:, [150]]=M_imp(beta[:,[150]], area[:,[150]])-M_evap(Twall[:,[150]],h_air[:,[150]],pressure[:,[150]], area[:,[150]]) #计算流出质量
	for i in range(150, 308):  #翼型上部
		Min[:,[i+1]]=Mout[:,[i]]   #水膜流动关系
		Mout[:,[i+1]]=M_imp(beta[:,[i+1]], area[:,[i]])+Min[:,[i+1]]-M_evap(Twall[:,[i+1]],h_air[:,[i+1]],pressure[:,[i+1]], area[:,[i+1]])
		#防冰状态下翼型上部流动质量关系
		if Mout[:,[i+1]]<0:#限制流出质量非负
			Mout[:,[i+1]]=0

	for i in range(0, 151):#翼型下部计算结果
		Min[:,[150-i-1]]=Mout[:,[150-i]]
		Mout[:,[150-i-1]]=M_imp(beta[:,[150-i-1]], area[:,[i]])+Min[:,[150-i-1]]-M_evap(Twall[:,[150-i-1]],h_air[:,[150-i-1]],pressure[:,[150-i-1]], area[:,[150-i-1]])
		if Mout[:,[150-i-1]]<0:#限制流出质量非负
			Mout[:,[150-i-1]]=0
	return Mout, Min

def temperature_i0(T):#the temperature after iteration of the first cell(leading edge)
	Mimp=M_imp(beta[0, 151], area[0, 151])  #impingment water mass
	Qcond=Q_cond(1586614.173, 1708661.417 , 2562992.126, 3417322.835, 2074803.15, 1464566.929, 1464566.929)#heat conduction
	Qaero=Q_aero(h_air[0, 151])*area[0, 151]#aerodynamic heat
	Qimp=Q_imp(beta[0, 151], area[0, 151])#energy in impinged water
	Mevap=M_evap(T, h_air[0, 151], pressure[0, 151], area[0, 151])
	Mout = Mimp-Mevap
	print(Mevap, Mout, Mimp)
	Qevap=Q_evap(Mevap, T, area[0,151])
	print(Qevap/area[0, 151])
	Twall_local=(Qcond[0, 151] * area[0, 151] +Qimp+Qaero-Qevap + h_air[0, 151]* T_inlet*area[0, 151])/(Mout * Cp_water + h_air[0, 151] * area[0, 151])
	return Twall_local
Twall0=Twall[0, 151]
print(temperature_i0(Twall0))






