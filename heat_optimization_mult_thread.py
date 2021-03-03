#In this code, I try to modified the single process to multi-process, but it is failed.
import matplotlib.pyplot as plt
import numpy as np
import geatpy as ea
import multiprocessing as mp
from multiprocessing import Pool as ProcessPool
from multiprocessing.dummy import Pool as ThreadPool




Cp_water=4200
V_inlet=89.4
T_inlet=251.35
Cp_air=1004.5
LWC=0.55
NIND=100
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
Twall=np.ones((1,309))*309#转换成矩阵
# print(Twall)

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
	def __init__(self, PoolType):
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
		#设置多线程还是多进程
		self.PoolType=PoolType
		if self.PoolType =='Thread':
			self.pool=ThreadPool(2) #设置池的大小
		elif self.PoolType == 'Process':
			num_cores=int(mp.cpu_count()) #统计计算机核心数
			self.pool=ProcessPool(num_cores)  #设置池的大小






	def aimFunc(self, pop):#目标函数


		Vars=pop.Phen#种群表现型，决策变量矩阵，每一列对应一个决策变量，每一行对应一个个体。




		for i in range(9):
			pop.ObjV=ProcessPool(9).map(aim2, Vars)
		#采用可行性法则处理约束，生成种群个体违反约束程度矩阵
		#违反约束程度矩阵，每一行对应种群的一个个体，每一列对应一个约束条件。
		#np.hstack把列向量拼合在一起形成约束矩阵。
		pop.CV=np.hstack([
			              # M_imp(beta[:, 306 ], area[0, 306]) + Mout[0, 305] - M_evap(x[:, [306]], h_air[0, 306], pressure[0, 306],  area[0, 306]),
						  # M_imp(beta[0, 307 ], area[0, 307]) + Mout[0, 306] - M_evap(x[:, [307]], h_air[0, 307],pressure[0, 307],	area[0, 307]),
						  	-(M_imp(beta[:, [308]], area[:, [308]]) + Mout[:, [307]] - M_evap(x[:, [308]], h_air[:, [308]], pressure[:, [308]],area[:, [308]])),
							M_imp(beta[:, [308] ], area[:, [308]]) + Mout[:, [307]] - M_evap(x[:, [308]], h_air[:, [308]],pressure[:, [308]],	area[:, [308]])-0.00002,
							-(M_imp(beta[:, [307]], area[:, [307]]) + Mout[:, [306]] - M_evap(x[:, [307]], h_air[:, [307]], pressure[:, [307]],area[:, [307]])),
							M_imp(beta[:, [307] ], area[:, [307]]) + Mout[:, [306]] - M_evap(x[:, [307]], h_air[:, [307]],pressure[:, [307]],	area[:, [307]])-0.00002,
							-(M_imp(beta[:, [306]], area[:, [306]]) + Mout[:, [305]] - M_evap(x[:, [306]], h_air[:, [306]], pressure[:, [306]],area[:, [306]])),
							M_imp(beta[:, [306] ], area[:, [306]]) + Mout[:, [305]] - M_evap(x[:, [306]], h_air[:, [306]],pressure[:, [306]],	area[:, [306]])-0.00002,
							-(M_imp(beta[:, [305]], area[:, [305]]) + Mout[:, [304]] - M_evap(x[:, [305]], h_air[:, [305]], pressure[:, [305]],area[:, [305]])),
							M_imp(beta[:, [305] ], area[:, [305]]) + Mout[:, [304]] - M_evap(x[:, [305]], h_air[:, [305]],pressure[:, [305]],	area[:, [305]])-0.00002,
							-(M_imp(beta[:, [304]], area[:, [304]]) + Mout[:, [303]] - M_evap(x[:, [304]], h_air[:, [304]], pressure[:, [304]],area[:, [304]])),
							M_imp(beta[:, [304] ], area[:, [304]]) + Mout[:, [303]] - M_evap(x[:, [304]], h_air[:, [304]],pressure[:, [304]],	area[:, [304]])-0.00002,
							-(M_imp(beta[:, [303]], area[:, [303]]) + Mout[:, [302]] - M_evap(x[:, [303]], h_air[:, [303]], pressure[:, [303]],area[:, [303]])),
							M_imp(beta[:, [303] ], area[:, [303]]) + Mout[:, [302]] - M_evap(x[:, [303]], h_air[:, [303]],pressure[:, [303]],	area[:, [303]])-0.00002,
							-(M_imp(beta[:, [302]], area[:, [302]]) + Mout[:, [301]] - M_evap(x[:, [302]], h_air[:, [302]], pressure[:, [302]],area[:, [302]])),
							M_imp(beta[:, [302] ], area[:, [302]]) + Mout[:, [301]] - M_evap(x[:, [302]], h_air[:, [302]],pressure[:, [302]],	area[:, [302]])-0.00002,
							-(M_imp(beta[:, [301]], area[:, [301]]) + Mout[:, [300]] - M_evap(x[:, [301]], h_air[:, [301]], pressure[:, [301]],area[:, [301]])),
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
							# -(M_imp(beta[:, [292]], area[:, [292]]) + Mout[:, [291]] - M_evap(x[:, [292]], h_air[:, [292]], pressure[:, [292]],area[:, [292]])),
							# M_imp(beta[:, [292] ], area[:, [292]]) + Mout[:, [291]] - M_evap(x[:, [292]], h_air[:, [292]],pressure[:, [292]],	area[:, [292]])-0.00002,
							# -(M_imp(beta[:, [291]], area[:, [291]]) + Mout[:, [290]] - M_evap(x[:, [291]], h_air[:, [291]], pressure[:, [291]],area[:, [291]])),
							# M_imp(beta[:, [291] ], area[:, [291]]) + Mout[:, [290]] - M_evap(x[:, [291]], h_air[:, [291]],pressure[:, [291]],	area[:, [291]])-0.00002,
							# -(M_imp(beta[:, [290]], area[:, [290]]) + Mout[:, [289]] - M_evap(x[:, [290]], h_air[:, [290]], pressure[:, [290]],area[:, [290]])),
							# M_imp(beta[:, [290] ], area[:, [290]]) + Mout[:, [289]] - M_evap(x[:, [290]], h_air[:, [290]],pressure[:, [290]],	area[:, [290]])-0.00002,
							# -(M_imp(beta[:, [289]], area[:, [289]]) + Mout[:, [288]] - M_evap(x[:, [289]], h_air[:, [289]], pressure[:, [289]],area[:, [289]])),
							# M_imp(beta[:, [289] ], area[:, [289]]) + Mout[:, [288]] - M_evap(x[:, [289]], h_air[:, [289]],pressure[:, [289]],	area[:, [289]])-0.00002,
			              #第一个约束，第一个网格内无runback water
						  -(M_imp(beta[:, [0]], area[:, [0]]) + Mout[:, [1]] - M_evap(x[:, [0]], h_air[:, [0]], pressure[:, [0]],   area[:, [0]])),#第二个约束，最后一个网格内无runback water
			              M_imp(beta[:, [0]], area[:, [0]]) + Mout[:, [1]] - M_evap(x[:, [0]], h_air[:, [0]], pressure[:, [0]], area[:, [0]])-0.000002,
						  np.abs(x[:,[307]]-x[:,[308]])-2,#第三个约束，最后一点和倒数第二点的温度差绝对值小于5
						  x[:, [150]]-310,
						  ])


		# print( M_evap(x[:, [308]], h_air[0, 308],pressure[0, 308],area[0, 308]))
		# print( M_imp(beta[0, 308 ], area[0, 308]))
		# print(x[-1, 308])
		# print(M_imp(beta[0, 150], area[0, 150]))
		# print(M_evap(x[-1, 150], h_air[0, 150], pressure[0, 150],  area[0, 150]))
		# y=np.arange(0, 309).reshape(1, 309)
		# plt.plot(y, Mout[[-1],:])
		# plt.show()

			def aim2(args):
				aim2 = 0
				heat_flux_array = np.zeros((NIND, 309))
				Mevap = np.zeros((NIND, 309))  # 蒸发质量初始化
				Qevap = np.zeros((NIND, 309))  # 蒸发能量初始化
				x = np.zeros((NIND, 309))  # 309cells数， NIND种群规模
				Min = np.zeros((NIND, 309))  # 流入质量
				Mout = np.zeros((NIND, 309))  # 流入质量
				for i in range(0, 309):
					x[:, [i]] = args[:, [i]]  # 每一行是一条染色体，也就是一个个体，一个具体问题的解
				for j in range(0, 309):
					Mevap[:, [j]] = M_evap(x[:, [j]], h_air[:, [j]], pressure[:, [j]],
										   area[:, [j]])  # 根据决策变量（也就是温度）计算蒸发质量
					Qevap[:, [j]] = Q_evap(M_evap(x[:, [j]], h_air[:, [j]], pressure[:, [j]], area[:, [j]]), x[:, [j]],
										   area[:, [j]])  # 计算蒸发能量
				# print(x[-1,[307]])
				Mout[:, [150]] = M_imp(beta[:, [150]], area[:, [150]]) - M_evap(x[:, [150]], h_air[:, [150]],
																				pressure[:, [150]],
																				area[:, [150]])  # 计算中间网格的流出量
				# print(Mout[0, 150])
				# print(x[-1, 150])
				for i in range(150, 308):  # mout和min没有进行进化，只是利用了种群中的一个个体，这是有问题滴
					Min[:, [i + 1]] = Mout[:, [i]]  # 在当前个体(温度的解)下获得水膜流动关系(翼型上部)
					Mout[:, [i + 1]] = M_imp(beta[:, [i + 1]], area[:, [i + 1]]) + Min[:, [i + 1]] - M_evap(
						x[:, [i + 1]], h_air[:, [i + 1]], pressure[:, [i + 1]], area[:, [i + 1]])
					if Mout[:, [i + 1]].all() < 0:  # 限制流出质量非负
						Mout[:, [i + 1]] = 0

				for i in range(0, 151):
					Min[:, [150 - i - 1]] = Mout[:, [150 - i]]  # 在当前个体(温度的解)下获得水膜流动关系(翼型下部)
					Mout[:, [150 - i - 1]] = M_imp(beta[:, [150 - i - 1]], area[:, [i]]) + Min[:,
																						   [150 - i - 1]] - M_evap(
						x[:, [150 - i - 1]], h_air[:, [150 - i - 1]], pressure[:, [150 - i - 1]],
						area[:, [150 - i - 1]])
					if Mout[:, [150 - i - 1]].all() < 0:  # 限制流出质量非负
						Mout[:, [150 - i - 1]] = 0

				for i in range(0, 309):
					heat_flux_array[:, [i]] = aim(x[:, [i]], Mout[:, [i]], Min[:, [i]], area[:, [i]], h_air[:, [i]],
												  pressure[:, [i]], beta[:, [i]]) * area[:, [i]]
					aim2 += heat_flux_array[:, [i]]
				return aim2
