import matplotlib.pyplot as plt
import numpy as np
import geatpy as ea
import math as m

#----------空气物性参数1--------------

R_a           = 287
Cp_air        = 1013

#----------------试验参数------------------
LWC           = 0.55*1e-3
MVD           = 20*1e-6
p_tot         = 70000
u_inf         = 89.4
T_tot         = 251.55
T0            = 273.15
T_inf         = T_tot-u_inf*u_inf/2/Cp_air
delt_s        = 1e-4
p_0           = p_tot/(1+(u_inf**2)/2/R_a/T_inf)


#--------------空气五行参数2----------------
rou_air       = p_0/R_a/T_inf
miu_air       = 1/(0.12764+124.38*(1/T_inf))*10**(-5)
lamda_air     = -1.4758e-2+2.3597e-3*(m.sqrt(T_inf))
Pr_air        = Cp_air*miu_air/lamda_air
M_air         = 29
R_air         = 287


#----------------液态水物性参数----------------------
lamda_water   = 0.599
M_H2O         = 18
R_v           = 461.4

Cp_water=4200
V_inlet=89.4
T_inlet=251.35
NIND=100
#I_water=2500000


#--------《Microphysics of Clouds and Precipitation》Magnus equation-----------
e_inf=610.70*m.exp(17.15*(T_inf-T0)/(T_inf-38.25)) #液膜表面饱和水蒸气压





#-----------读取参数--------------------
arealist=[]  #面积
betalist=[]   #水收集系数
h_airlist=[]   #空气对流换热系数
pressurelist=[]   #压力
Twall_steadylist=[]   #稳态壁温

#--------------------------------------
#-----------------读文件----------------
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
		for i in range(4, len(datasetlist)):#从第四行开始读
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
Twall=np.ones((1,309))*309#转换成矩阵
# print(Twall)








#------------------------------------------

def Q_extra(Twall, hair):#计算对流换热损失
	T_inlet=251.35
	Qextra=hair*(Twall-T_inlet)
	return Qextra

def Q_out(Twall, Mout, area):#计算流出热流密度
	Qout=Mout*Twall*Cp_water/area
	return Qout

def Q_in(Twall, Min, area):#计算流入热流密度
	Qin=Min*Twall*Cp_water/area
	return Qin

def Q_imp(beta,area):#计算撞击热流密度
	Qimp=(V_inlet * beta * 0.55 * 0.001*area)*(V_inlet**2)/2+(V_inlet * beta * LWC * 0.001*area)*Cp_water*T_inlet
	Qimp=Qimp/area
	return Qimp

def M_imp(beta, area):#计算撞击质量
	Mimp=(V_inlet * beta * 0.55 * 0.001*area)
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
	Qevap=(Mevap*I_water+Mevap*Twall*Cp_water)/area
	return Qevap


def Q_cond(Qout, Qextra, Qevap, Qimp, Qaero, Qin):#计算防冰负荷
	Qcond=Qout+Qextra+Qevap-Qimp-Qaero-Qin
	return Qcond




def Ite_water_film():
	T_s = np.zeros((1, 309))
	T_wall = np.zeros((1, 309))
	T_rec = np.zeros((1, 309))
	e_sat_water=np.zeros((1, 309))
	m_evap_per=np.zeros((1, 309))
	Le=np.zeros((1, 309))
	m_evap=np.zeros((1, 309))
	m_in=np.zeros((1, 309))
	for i in range(0, 151):
		if i==0:
			T_s[0, i]=273.15
			T_wall[0, i]=T_s[0, i]-T0
			for j in range(0, 999):
				T_ref=(T_s[0, i]+T_wall[0, i]+T0)
				if T_ref >= T0:
					row_water=(999.8396+18.224944*(T_ref-T0)-7.922210e-3*(T_ref-T0)**2.-55.44846e-6*(T_ref-T0)**3+149.7562e-9*(T_ref-T0)**4-393.2952e-12*(T_ref-T0)**5)/(1+(18.159725e-3)*(T_ref-T0))
					Cp_water = 4185.8518 * (0.9979 + 3.1e-6 * (T_ref - T0 - 35) ** 2 + 3.8e-9 * (T_ref - T0 - 35) ** 4)
					miu_water = 1e-3 * (1.76 * m.exp(-3.5254e-2 * (T_ref - T0) + 4.7163e-4 * (T_ref - T0) ** 2. - 6.0667e-6 * (T_ref - T0) ** 3))
				else:
					rou_water = 1000 * (0.99986 + 6.69e-5 * (T_ref - T0) - 8.486e-6 * (T_ref - T0) ** 2. + 1.518e-7 * (T_ref - T0) ** 3 - 6.9484e-9 * (T_ref - T0) ** 4 - 3.6449e-10 * (T_ref - T0) ** 5. - 7.497e-12 * (T_ref - T0) ** 6)
					Cp_water = 4185.8518 * (1.000938 - 2.7052e-3 * (T_ref - T0) - 2.3235e-5 * (T_ref - T0) ** 2. + 4.3778e-6 * (T_ref - T0) ** 3 + 2.7136e-7 * (T_ref - T0) ** 4)
					miu_water = 1e-3 * (1.76 * m.exp(-5.5721e-2 * (T_ref - T0) - 1.3943e-3 * (T_ref - T0) ** 2. - 4.3015e-5 * (T_ref - T0) ** 3))
				D_v= 0.211*((0.5*(T_ref+T_inf)/273.15)**1.94)*(101325/p_0)*1e-4
				Sc=miu_air/rou_air/D_v
				e_sat_water[0, i] = 610.70 * m.exp(17.15 * (T_s[0, i] - T0) / (T_s[0, i]- 38.25))

				Le[0, i] = 4185.8518 * (597.3 - 0.561 * (T_s[0, i] - T0))
				m_evap_per[0, i] = 0.622 * h_air[0, i] / Cp_air * ( e_sat_water[0, i] / (pressure[0, i] - e_sat_water[0, i]) - e_inf / (p_0 - e_inf))
				m_evap[0, i] = m_evap_per[0, i] * area[0, i]
				#-------------质量平衡-------------------------
				m_in[0, i]=0
				if (m_in[0, 1]+M_imp(beta[0,i], area[0, i])<=m_evap[0, i]):
					m_sta_p=0
					m_evap[0, i]=m_in[0, 1]+M_imp(beta[0,i], area[0, i])
					m_evap[0, i]=m_evap_per[0, i] *area[0, i]
































def Water_film():  #计算水膜
	Mevap=np.zeros((NIND, 309))   #初始化
	Qevap=np.zeros((NIND, 309))    #初始化
	for i in range(0, 309):
		Mevap[:,[i]]=M_evap(Twall[:,[i]], h_air[:,[i]], pressure[:,[i]], area[:,[i]])   #计算蒸发质量
		Qevap[:,i]=Q_evap(M_evap(Twall[:,[i]], h_air[:,[i]], pressure[:,[i]], area[:,[i]]), Twall[:,[i]], area[:,[i]])   #计算蒸发能量
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


def aim(Twall,Mout, Min, area, h_air, pressure, beta ):  # 目标函数
	return Q_out(Twall, Mout, area)+Q_extra(Twall, h_air)+Q_evap(M_evap(Twall, h_air, pressure, area),Twall, area)-Q_imp(beta, area)-Q_aero(h_air)-Q_in(Twall, Min, area)


class MyProblem(ea.Problem):#问题类
	def __init__(self):
		name='MyProblem'  #初始化函数名称
		M=1  #初始化目标维数
		maxormins=[1]  #初始化最大最小值，1：最小化该目标；-1：最大化该目标
		Dim=309 #初始化决策变量维数
		varTypes=[0]*Dim #初始化决策变量类型，元素为0表示对应的变量是连续的；1表示是离散的，连续型变量
		lb=[290 for x in range(0, 309)] #决策变量下界
		ub=[320 for x in range(0, 309)]  #决策变量上界
		lbin=[0 for x in range(0,309)]  ## 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
		ubin=[1 for x in range(0, 309)]  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
		#调用父类构造方法完成实例化
		ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
	def aimFunc(self, pop):#目标函数
		aim2=0
		heat_flux_array_area = np.zeros((NIND, 309))
		heat_flux_array = np.zeros((NIND, 309))
		Mevap = np.zeros((NIND, 309)) #蒸发质量初始化
		Qevap = np.zeros((NIND, 309)) #蒸发能量初始化

		Vars=pop.Phen#种群表现型，决策变量矩阵，每一列对应一个决策变量，每一行对应一个个体。
		x = np.zeros((NIND, 309))#309cells数， NIND种群规模
		Min=np.zeros((NIND,309))#流入质量
		Mout=np.zeros((NIND,309))#流入质量
		for i in range(0,309):
			x[:,[i]] = Vars[:, [i]] #每一行是一条染色体，也就是一个个体，一个具体问题的解
		for j in range(0, 309):
			Mevap[:, [j]] = M_evap(x[:,[j]], h_air[:, [j]], pressure[:, [j]], area[:, [j]])#根据决策变量（也就是温度）计算蒸发质量
			Qevap[:, [j]] = Q_evap(M_evap(x[:,[j]],  h_air[:, [j]], pressure[:, [j]], area[:, [j]]), x[:,[j]], area[:, [j]])#计算蒸发能量
		# print(x[-1,[307]])
		Mout[:, [150]] = M_imp(beta[:, [150]], area[:, [150]]) - M_evap(x[:, [150]], h_air[:, [150]], pressure[:, [150]], area[:, [150]])#计算中间网格的流出量
		# print(Mout[0, 150])
		# print(x[-1, 150])
		for i in range(150, 308):#mout和min没有进行进化，只是利用了种群中的一个个体，这是有问题滴
			Min[:, [i + 1]] = Mout[:, [i]]#在当前个体(温度的解)下获得水膜流动关系(翼型上部)
			Mout[:, [i + 1]] = M_imp(beta[:, [i + 1]], area[:, [i+1]]) + Min[:,[ i + 1]] - M_evap(x[:, [i+1]], h_air[:, [i+1]],pressure[:, [i+1]],	area[:, [i+1]])
			if Mout[:, [i + 1]].all() < 0:#限制流出质量非负
				Mout[:, [i + 1]] = 0

		for i in range(0, 151):
			Min[:, [150 - i - 1]] = Mout[:, [150 - i]]#在当前个体(温度的解)下获得水膜流动关系(翼型下部)
			Mout[:, [150 - i - 1]] = M_imp(beta[:, [150 - i - 1]], area[:, [i]]) + Min[:, [150 - i - 1]] - M_evap(x[:,[150-i-1]], h_air[:, [150 - i - 1]], pressure[:, [150 - i - 1]], area[:, [150 - i - 1]])
			if Mout[:, [150 - i - 1]].all() < 0:#限制流出质量非负
				Mout[:, [150 - i - 1]] = 0



		for i in range(0, 309):
			heat_flux_array_area[:,[i]]=aim(x[:, [i]], Mout[:,[i]], Min[:,[i]], area[:,[i]],h_air[:,[i]], pressure[:,[i]], beta[:,[i]])*area[:, [i]]
			heat_flux_array[:,[i]]=aim(x[:, [i]], Mout[:,[i]], Min[:,[i]], area[:,[i]],h_air[:,[i]], pressure[:,[i]], beta[:,[i]])
			aim2+=heat_flux_array_area[:,[i]]
		pop.ObjV=aim2
		#采用可行性法则处理约束，生成种群个体违反约束程度矩阵
		#违反约束程度矩阵，每一行对应种群的一个个体，每一列对应一个约束条件。
		#np.hstack把列向量拼合在一起形成约束矩阵。
		pop.CV=np.hstack([
			              # M_imp(beta[:, 306 ], area[0, 306]) + Mout[0, 305] - M_evap(x[:, [306]], h_air[0, 306], pressure[0, 306],  area[0, 306]),
						  # M_imp(beta[0, 307 ], area[0, 307]) + Mout[0, 306] - M_evap(x[:, [307]], h_air[0, 307],pressure[0, 307],	area[0, 307]),
						  -(M_imp(beta[:, [308]], area[:, [308]]) + Mout[:, [307]] - M_evap(x[:, [308]], h_air[:, [308]], pressure[:, [308]],area[:, [308]])),
							M_imp(beta[:, [308] ], area[:, [308]]) + Mout[:, [307]] - M_evap(x[:, [308]], h_air[:, [308]],pressure[:, [308]],	area[:, [308]])-0.00002,
							# -(M_imp(beta[:, [307]], area[:, [307]]) + Mout[:, [306]] - M_evap(x[:, [307]], h_air[:, [307]], pressure[:, [307]],area[:, [307]])),
							# M_imp(beta[:, [307] ], area[:, [307]]) + Mout[:, [306]] - M_evap(x[:, [307]], h_air[:, [307]],pressure[:, [307]],	area[:, [307]])-0.00002,
							# -(M_imp(beta[:, [306]], area[:, [306]]) + Mout[:, [305]] - M_evap(x[:, [306]], h_air[:, [306]], pressure[:, [306]],area[:, [306]])),
							# M_imp(beta[:, [306] ], area[:, [306]]) + Mout[:, [305]] - M_evap(x[:, [306]], h_air[:, [306]],pressure[:, [306]],	area[:, [306]])-0.00002,
							# -(M_imp(beta[:, [305]], area[:, [305]]) + Mout[:, [304]] - M_evap(x[:, [305]], h_air[:, [305]], pressure[:, [305]],area[:, [305]])),
							# M_imp(beta[:, [305] ], area[:, [305]]) + Mout[:, [304]] - M_evap(x[:, [305]], h_air[:, [305]],pressure[:, [305]],	area[:, [305]])-0.00002,
							# -(M_imp(beta[:, [304]], area[:, [304]]) + Mout[:, [303]] - M_evap(x[:, [304]], h_air[:, [304]], pressure[:, [304]],area[:, [304]])),
							# M_imp(beta[:, [304] ], area[:, [304]]) + Mout[:, [303]] - M_evap(x[:, [304]], h_air[:, [304]],pressure[:, [304]],	area[:, [304]])-0.00002,
							# -(M_imp(beta[:, [303]], area[:, [303]]) + Mout[:, [302]] - M_evap(x[:, [303]], h_air[:, [303]], pressure[:, [303]],area[:, [303]])),
							# M_imp(beta[:, [303] ], area[:, [303]]) + Mout[:, [302]] - M_evap(x[:, [303]], h_air[:, [303]],pressure[:, [303]],	area[:, [303]])-0.00002,
							# -(M_imp(beta[:, [302]], area[:, [302]]) + Mout[:, [301]] - M_evap(x[:, [302]], h_air[:, [302]], pressure[:, [302]],area[:, [302]])),
							# M_imp(beta[:, [302] ], area[:, [302]]) + Mout[:, [301]] - M_evap(x[:, [302]], h_air[:, [302]],pressure[:, [302]],	area[:, [302]])-0.00002,
							# -(M_imp(beta[:, [301]], area[:, [301]]) + Mout[:, [300]] - M_evap(x[:, [301]], h_air[:, [301]], pressure[:, [301]],area[:, [301]])),
							# M_imp(beta[:, [301] ], area[:, [301]]) + Mout[:, [300]] - M_evap(x[:, [301]], h_air[:, [301]],pressure[:, [301]],	area[:, [301]])-0.00002,
							# -(M_imp(beta[:, [300]], area[:, [300]]) + Mout[:, [299]] - M_evap(x[:, [300]], h_air[:, [300]], pressure[:, [300]],area[:, [300]])),
							# M_imp(beta[:, [300] ], area[:, [300]]) + Mout[:, [299]] - M_evap(x[:, [300]], h_air[:, [300]],pressure[:, [300]],	area[:, [300]])-0.00002,
							# -(M_imp(beta[:, [299]], area[:, [299]]) + Mout[:, [298]] - M_evap(x[:, [299]], h_air[:, [299]], pressure[:, [299]],area[:, [299]])),
							# M_imp(beta[:, [299] ], area[:, [299]]) + Mout[:, [298]] - M_evap(x[:, [299]], h_air[:, [299]],pressure[:, [299]],	area[:, [299]])-0.00002,
							# -(M_imp(beta[:, [298]], area[:, [298]]) + Mout[:, [297]] - M_evap(x[:, [298]], h_air[:, [298]], pressure[:, [298]],area[:, [298]])),
							# M_imp(beta[:, [298] ], area[:, [298]]) + Mout[:, [297]] - M_evap(x[:, [298]], h_air[:, [298]],pressure[:, [298]],	area[:, [298]])-0.00002,
							# -(M_imp(beta[:, [297]], area[:, [297]]) + Mout[:, [296]] - M_evap(x[:, [297]], h_air[:, [297]], pressure[:, [297]],area[:, [297]])),
							# M_imp(beta[:, [297] ], area[:, [297]]) + Mout[:, [296]] - M_evap(x[:, [297]], h_air[:, [297]],pressure[:, [297]],	area[:, [297]])-0.00002,
							# -(M_imp(beta[:, [296]], area[:, [296]]) + Mout[:, [295]] - M_evap(x[:, [296]], h_air[:, [296]], pressure[:, [296]],area[:, [296]])),
							# M_imp(beta[:, [296] ], area[:, [296]]) + Mout[:, [295]] - M_evap(x[:, [296]], h_air[:, [296]],pressure[:, [296]],	area[:, [296]])-0.00002,
							# -(M_imp(beta[:, [295]], area[:, [295]]) + Mout[:, [294]] - M_evap(x[:, [295]], h_air[:, [295]], pressure[:, [295]],area[:, [295]])),
							# M_imp(beta[:, [295] ], area[:, [295]]) + Mout[:, [294]] - M_evap(x[:, [295]], h_air[:, [295]],pressure[:, [295]],	area[:, [295]])-0.00002,
							# -(M_imp(beta[:, [294]], area[:, [294]]) + Mout[:, [293]] - M_evap(x[:, [294]], h_air[:, [294]], pressure[:, [294]],area[:, [294]])),
							# M_imp(beta[:, [294] ], area[:, [294]]) + Mout[:, [293]] - M_evap(x[:, [294]], h_air[:, [294]],pressure[:, [294]],	area[:, [294]])-0.00002,
							# -(M_imp(beta[:, [293]], area[:, [293]]) + Mout[:, [292]] - M_evap(x[:, [293]], h_air[:, [293]], pressure[:, [293]],area[:, [293]])),
							# M_imp(beta[:, [293] ], area[:, [293]]) + Mout[:, [292]] - M_evap(x[:, [293]], h_air[:, [293]],pressure[:, [293]],	area[:, [293]])-0.00002,
							#-(M_imp(beta[:, [292]], area[:, [292]]) + Mout[:, [291]] - M_evap(x[:, [292]], h_air[:, [292]], pressure[:, [292]],area[:, [292]])),
							# M_imp(beta[:, [292] ], area[:, [292]]) + Mout[:, [291]] - M_evap(x[:, [292]], h_air[:, [292]],pressure[:, [292]],	area[:, [292]])-0.00002,
							# -(M_imp(beta[:, [291]], area[:, [291]]) + Mout[:, [290]] - M_evap(x[:, [291]], h_air[:, [291]], pressure[:, [291]],area[:, [291]])),
							# M_imp(beta[:, [291] ], area[:, [291]]) + Mout[:, [290]] - M_evap(x[:, [291]], h_air[:, [291]],pressure[:, [291]],	area[:, [291]])-0.00002,
							# -(M_imp(beta[:, [290]], area[:, [290]]) + Mout[:, [289]] - M_evap(x[:, [290]], h_air[:, [290]], pressure[:, [290]],area[:, [290]])),
							# M_imp(beta[:, [290] ], area[:, [290]]) + Mout[:, [289]] - M_evap(x[:, [290]], h_air[:, [290]],pressure[:, [290]],	area[:, [290]])-0.00002,
							# -(M_imp(beta[:, [289]], area[:, [289]]) + Mout[:, [288]] - M_evap(x[:, [289]], h_air[:, [289]], pressure[:, [289]],area[:, [289]])),
							# M_imp(beta[:, [289] ], area[:, [289]]) + Mout[:, [288]] - M_evap(x[:, [289]], h_air[:, [289]],pressure[:, [289]],	area[:, [289]])-0.00002,
			              #第一个约束，第一个网格内无runback water
						 	-(M_imp(beta[:, [0]], area[:, [0]]) + Mout[:, [1]] - M_evap(x[:, [0]], h_air[:, [0]], pressure[:, [0]],   area[:, [0]])),#最后一个网格溢流水大于0
			             	M_imp(beta[:, [0]], area[:, [0]]) + Mout[:, [1]] - M_evap(x[:, [0]], h_air[:, [0]], pressure[:, [0]], area[:, [0]])-0.000002,#最后一个网格溢流水小于0.000002
						   	np.abs(x[:,[307]]-x[:,[308]])-2,#第三个约束，最后一点和倒数第二点的温度差绝对值小于2
						  # x[:, [150]]-310,
			              	np.abs(heat_flux_array[:, [0]]-heat_flux_array[:, [1]])-1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [2]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [3]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [4]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [5]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [6]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [7]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [8]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [9]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [10]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [11]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [12]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [13]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [14]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [15]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [16]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [17]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [18]]) - 1000,
							np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [19]]) - 1000,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [21]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [22]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [23]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [24]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [25]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [26]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [27]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [28]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [29]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [30]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [31]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [32]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [33]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [34]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [35]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [36]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [37]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [38]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [39]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [40]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [41]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [42]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [43]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [44]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [45]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [46]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [47]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [48]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [49]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [50]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [51]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [52]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [53]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [54]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [55]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [56]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [57]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [58]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [59]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [60]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [61]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [62]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [63]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [64]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [65]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [66]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [67]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [68]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [69]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [70]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [71]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [72]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [73]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [74]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [75]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [76]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [77]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [78]]) - 10,
							# np.abs(heat_flux_array[:, [0]] - heat_flux_array[:, [79]]) - 10,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [101]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [102]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [103]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [104]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [105]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [106]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [107]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [108]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [109]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [110]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [111]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [112]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [113]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [114]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [115]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [116]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [117]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [118]]) - 1000,
							np.abs(heat_flux_array[:, [100]] - heat_flux_array[:, [119]]) - 1000,
							np.abs(heat_flux_array[:, [200]] - heat_flux_array[:, [201]]) - 1000,
							np.abs(heat_flux_array[:, [200]] - heat_flux_array[:, [202]]) - 1000,
							np.abs(heat_flux_array[:, [200]] - heat_flux_array[:, [203]]) - 1000,
							np.abs(heat_flux_array[:, [200]] - heat_flux_array[:, [204]]) - 1000,
							np.abs(heat_flux_array[:, [200]] - heat_flux_array[:, [205]]) - 1000,
							np.abs(heat_flux_array[:, [200]] - heat_flux_array[:, [206]]) - 1000,
							np.abs(heat_flux_array[:, [200]] - heat_flux_array[:, [207]]) - 1000,
							np.abs(heat_flux_array[:, [200]] - heat_flux_array[:, [208]]) - 1000,
							np.abs(heat_flux_array[:, [200]] - heat_flux_array[:, [209]]) - 1000,
		                  ])
		# print( M_evap(x[:, [308]], h_air[0, 308],pressure[0, 308],area[0, 308]))
		# print( M_imp(beta[0, 308 ], area[0, 308]))
		# print(x[-1, 308])
		# print(M_imp(beta[0, 150], area[0, 150]))
		# print(M_evap(x[-1, 150], h_air[0, 150], pressure[0, 150],  area[0, 150]))
		# y=np.arange(0, 309).reshape(1, 309)
		# plt.plot(y, Mout[[-1],:])
		# plt.show()


