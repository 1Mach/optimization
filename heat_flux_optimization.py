import matplotlib.pyplot as plt
import numpy as np
Cp_water=4200
V_inlet=89.4
T_inlet=251.35
Cp_air=1004.5
LWC=0.55
I_water=2500000
area=[]
beta=[]
h_air=[]
pressure=[]
Twall_steady=[]
def read_original_data(filename):
	s = []
	area = []
	dataset = []
	dataset2 = []
	dataset3 = []
	with open(filename, 'r') as f:
		lines = f.readlines()
		for data in lines:
			lines1 = data.strip('\n')
			dataset.append(lines1)
		for i in range(4, len(dataset)):
			dataset2.append(dataset[i])
		for j in range(0, len(dataset2)):
			dataset3.append(dataset2[j].split('\t'))
		for i in range(0, len(dataset3) - 1):
			s.append(float(dataset3[i][0]))
			area.append(float(dataset3[i][1]))
	return area

area.append(read_original_data('./original_data/area'))
beta.append(read_original_data('./original_data/beta'))
h_air.append(read_original_data('./original_data/heat_transfer_coefficient'))
pressure.append(read_original_data('./original_data/pressure'))
Twall_steady.append(read_original_data('./original_data/temperature'))
beta = beta[0]
area = area[0]
h_air = h_air[0]
pressure=pressure[0]
Twall_steady=Twall_steady[0]
Twall1=[300 for x1 in range(0,309)]
Twall2=[306 for x2 in range(0,309)]
Twall3=[308 for x3 in range(0,309)]

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

Mevap=[]
Qevap=[]
for i in range(0, 309):
	Mevap.append(M_evap(Twall_steady[i], h_air[i], pressure[i], area[i]))
	Qevap.append(Q_evap(M_evap(Twall_steady[i], h_air[i], pressure[i], area[i]), Twall_steady[i], area[i]))
Mout=[0 for x in range(0,309)]
Min=[0 for x in range(0,309)]
Mout[150]=M_imp(beta[150], area[150])-M_evap(Twall_steady[150],h_air[150],pressure[150], area[150])
for i in range(150, 308):
	Min[i+1]=Mout[i]
	Mout[i+1]=M_imp(beta[i+1], area[i])+Min[i+1]-M_evap(Twall_steady[i+1],h_air[i+1],pressure[i+1], area[i+1])
	if Mout[i+1]<0:
		Mout[i+1]=0
for i in range(0, 151):
	Min[150-i-1]=Mout[150-i]
	Mout[150-i-1]=M_imp(beta[150-i-1], area[i])+Min[150-i-1]-M_evap(Twall_steady[150-i-1],h_air[150-i-1],pressure[150-i-1], area[150-i-1])
	if Mout[150-i-1]<0:
		Mout[150-i-1]=0

Mevap1=[]
Qevap1=[]
for i in range(0, 309):
	Mevap1.append(M_evap(Twall1[i], h_air[i], pressure[i], area[i]))
	Qevap1.append(Q_evap(M_evap(Twall1[i], h_air[i], pressure[i], area[i]), Twall1[i], area[i]))
Mout1=[0 for x in range(0,309)]
Min1=[0 for x in range(0,309)]
Mout1[150]=M_imp(beta[150], area[150])-M_evap(Twall1[150],h_air[150],pressure[150], area[150])
for i in range(150, 308):
	Min1[i+1]=Mout1[i]
	Mout1[i+1]=M_imp(beta[i+1], area[i])+Min1[i+1]-M_evap(Twall1[i+1],h_air[i+1],pressure[i+1], area[i+1])
	if Mout1[i+1]<0:
		Mout1[i+1]=0
for i in range(0, 151):
	Min1[150-i-1]=Mout1[150-i]
	Mout1[150-i-1]=M_imp(beta[150-i-1], area[i])+Min1[150-i-1]-M_evap(Twall1[150-i-1],h_air[150-i-1],pressure[150-i-1], area[150-i-1])
	if Mout1[150-i-1]<0:
		Mout1[150-i-1]=0

Mevap2=[]
Qevap2=[]
for i in range(0, 309):
	Mevap2.append(M_evap(Twall2[i], h_air[i], pressure[i], area[i]))
	Qevap2.append(Q_evap(M_evap(Twall2[i], h_air[i], pressure[i], area[i]), Twall2[i], area[i]))
Mout2=[0 for x in range(0,309)]
Min2=[0 for x in range(0,309)]
Mout2[150]=M_imp(beta[150], area[150])-M_evap(Twall2[150],h_air[150],pressure[150], area[150])
for i in range(150, 308):
	Min2[i+1]=Mout2[i]
	Mout2[i+1]=M_imp(beta[i+1], area[i])+Min2[i+1]-M_evap(Twall2[i+1],h_air[i+1],pressure[i+1], area[i+1])
	if Mout2[i+1]<0:
		Mout2[i+1]=0

for i in range(0, 151):
	Min2[150-i-1]=Mout2[150-i]
	Mout2[150-i-1]=M_imp(beta[150-i-1], area[i])+Min2[150-i-1]-M_evap(Twall2[150-i-1],h_air[150-i-1],pressure[150-i-1], area[150-i-1])
	if Mout2[150-i-1]<0:
		Mout2[150-i-1]=0

Mevap3=[]
Qevap3=[]
for i in range(0, 309):
	Mevap3.append(M_evap(Twall3[i], h_air[i], pressure[i], area[i]))
	Qevap3.append(Q_evap(M_evap(Twall3[i], h_air[i], pressure[i], area[i]), Twall3[i], area[i]))
Mout3=[0 for x in range(0,309)]
Min3=[0 for x in range(0,309)]
Mout3[150]=M_imp(beta[150], area[150])-M_evap(Twall3[150],h_air[150],pressure[150], area[150])
for i in range(150, 308):
	Min3[i+1]=Mout3[i]
	Mout3[i+1]=M_imp(beta[i+1], area[i])+Min3[i+1]-M_evap(Twall3[i+1],h_air[i+1],pressure[i+1], area[i+1])
	if Mout3[i+1]<0:
		Mout3[i+1]=0

for i in range(0, 151):
	Min3[150-i-1]=Mout3[150-i]
	Mout3[150-i-1]=M_imp(beta[150-i-1], area[i])+Min3[150-i-1]-M_evap(Twall3[150-i-1],h_air[150-i-1],pressure[150-i-1], area[150-i-1])
	if Mout3[150-i-1]<0:
		Mout3[150-i-1]=0
# def aim(Twall):  # 传入种群染色体矩阵解码后的基因表现型矩阵
# 	return Q_out(Twall, Mout[10], area[10])+Q_extra(Twall, h_air[10])+Q_evap(M_evap(Twall, h_air[10], pressure[10], area[10]),Twall, area[10])-Q_imp(beta[10], area[10])-Q_aero(h_air[10])-Q_in(Twall, Min[10], area[10])
Qcond=[]
for i in range(0, 309):
	Qcond_i=Q_out(Twall_steady[i], Mout[i], area[i])+Q_extra(Twall_steady[i], h_air[i])+Q_evap(M_evap(Twall_steady[i], h_air[i], pressure[i], area[i]),Twall_steady[i], area[i])-Q_imp(beta[i], area[i])-Q_aero(h_air[i])-Q_in(Twall_steady[i], Min[i], area[i])
	Qcond.append(Qcond_i)

Qcond1=[]
for i in range(0, 309):
	Qcond1_i=Q_out(Twall1[i], Mout1[i], area[i])+Q_extra(Twall1[i], h_air[i])+Q_evap(M_evap(Twall1[i], h_air[i], pressure[i], area[i]),Twall1[i], area[i])-Q_imp(beta[i], area[i])-Q_aero(h_air[i])-Q_in(Twall1[i], Min1[i], area[i])
	Qcond1.append(Qcond1_i)
Qcond2=[]
for i in range(0, 309):
	Qcond2_i=Q_out(Twall2[i], Mout2[i], area[i])+Q_extra(Twall2[i], h_air[i])+Q_evap(M_evap(Twall2[i], h_air[i], pressure[i], area[i]),Twall2[i], area[i])-Q_imp(beta[i], area[i])-Q_aero(h_air[i])-Q_in(Twall2[i], Min2[i], area[i])
	Qcond2.append(Qcond2_i)

Qcond3=[]
for i in range(0, 309):
	Qcond3_i=Q_out(Twall3[i], Mout3[i], area[i])+Q_extra(Twall3[i], h_air[i])+Q_evap(M_evap(Twall3[i], h_air[i], pressure[i], area[i]),Twall3[i], area[i])-Q_imp(beta[i], area[i])-Q_aero(h_air[i])-Q_in(Twall3[i], Min3[i], area[i])
	Qcond3.append(Qcond3_i)



x=np.arange(253.15, 273.15, 1)
y=np.arange(0, 309, 1)
y_list=y.tolist()
ax1=plt.subplot(211)
ax2=plt.subplot(212)
plt.sca(ax1)
# plt.plot(y, Min1)
# plt.xlabel('Number of cells')
plt.ylabel('Mout')
plt.plot(y, Mout,label='steady')
plt.plot(y, Mout1,label='300K')
plt.plot(y, Mout2,label='306K')
plt.plot(y, Mout3,label='308K')
plt.legend()
plt.sca(ax2)
plt.xlabel('Number of cells')
plt.ylabel('Qcond')
plt.plot(y, Qcond, label='steady')
plt.plot(y, Qcond1, label='300K')
plt.plot(y, Qcond2, label='306K')
plt.plot(y, Qcond3, label='308K')
# plt.plot(y, Twall_steady)
plt.legend()
plt.show()