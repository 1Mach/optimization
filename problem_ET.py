# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 19:22:39 2021

@author: guoxf
"""
import numpy as np
import geatpy as ea
import math as m

n=309
NIND=10
#-----------------空气参数----------------------
R_a           = 287
Cp_air        = 1013



#-----------------------实验参数---------------------------
LWC           = 0.55*1e-3
MVD           = 20*1e-6
p_tot         = 70000
u_inf         = 89.4
T_tot         = 251.55
T0            = np.ones((NIND, 1))*273.15
T_inf         = T_tot-u_inf*u_inf/2/Cp_air
delt_s        = 1e-4
p_0           = p_tot/(1+(u_inf**2)/2/R_a/T_inf)


#---------------------空气物性参数---------------------------
rou_air       = p_0/R_a/T_inf
miu_air       = 1/(0.12764+124.38*(1/T_inf))*10**(-5)
lamda_air     = -1.4758e-2+2.3597e-3*(m.sqrt(T_inf))
Pr_air        = Cp_air*miu_air/lamda_air
M_air         = 29
R_air         = 287



#--------------液态水物性-------------------------
lamda_water   = 0.599
M_H2O         = 18
R_v           = 461.4



#-----------------蒙皮导热系数-------------------
lamda_skin   = 4.7

#-----------------液膜表面饱和水蒸气压------------
e_inf=610.70*m.exp(17.15*(T_inf-273.15)/(T_inf-38.25)) #本来是(T_inf-T0(T0是冻结温度))



arealist=[]
betalist=[]
h_airlist=[]
pressurelist=[]
Twall_steadylist=[]
def read_original_data(filename):   #读取计算参数
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

betalist.append(read_original_data('./original_data/beta'))
beta=np.asarray(betalist)
beta=np.repeat(beta, NIND, axis=0)

h_airlist.append(read_original_data('./original_data/heat_transfer_coefficient'))
h_air=np.asarray(h_airlist)
h_air=np.repeat(h_air, NIND, axis=0)


pressurelist.append(read_original_data('./original_data/pressure'))
pressure=np.asarray(pressurelist)
pressure=np.repeat(pressure, NIND, axis=0)

T_wall=np.zeros((NIND, n))
T_s=np.zeros((NIND, n))
Le=np.zeros((NIND, n))
e_sat_water=np.zeros((NIND, n))
m_evap_per=np.zeros((NIND, n))
m_evap=np.zeros((NIND, n))
m_in=np.zeros((NIND, n))
m_out=np.zeros((NIND, n))
m_stat_p=np.zeros((NIND, n))
m_imp=np.zeros((NIND, n))
m_imp_per=np.zeros((NIND, n))
delt_w=np.zeros((NIND, n))
Q_in=np.zeros((NIND, n))
Q_out=np.zeros((NIND, n))
Q_imp=np.zeros((NIND, n))
Q_imp_per=np.zeros((NIND, n))
Q_evap=np.zeros((NIND, n))
Q_evap_per=np.zeros((NIND, n))
Q_anti=np.zeros((NIND, n))
Q_conv=np.zeros((NIND, n))
Q_conv_per=np.zeros((NIND, n))
Q_sta_p=np.zeros((NIND, n))
q_anti=np.zeros((NIND, n))

rou_water= np.zeros((NIND, 1))
Cp_water= np.zeros((NIND, 1))
miu_water=np.zeros((NIND, 1))
T_ref_1=np.zeros((NIND, 1))
T_ref_2=np.zeros((NIND, 1))
err=np.zeros((NIND, 1))
err_2=np.zeros((NIND, 1))


def cell_iteration(cell_num, NIND_num, T_wall_c, T_s_c):
	i_c = cell_num
	n_c = NIND_num
	T0_c=273.15

	T_ref_c = 0.5 * (T_s_c + T_wall_c + T0_c)
	if (T_ref_c >= T0_c):
		Cp_water = 4185.8518 * (0.9979 + 3.1e-6 * (T_ref_c - T0_c - 35) ** 2 + 3.8e-9 * (T_ref_c - T0_c - 35) ** 4)
	else:
		Cp_water = 4185.8518 * (1.000938 - 2.7052e-3 * (T_ref_c - T0_c) - 2.3235e-5 * (T_ref_c - T0_c) ** 2. + 4.3778e-6 * (T_ref_c - T0_c) ** 3 + 2.7136e-7 * (T_ref_c - T0_c) ** 4)
	e_sat_water_c = 610.70 * m.exp(17.15 * (T_s_c - T0_c) / (T_s_c - 38.25))
	Le_c = 4185.8518 * (597.3 - 0.561 * (T_s_c - T0_c))
	m_evap_per_c = 0.622 * h_air[0, i_c] / Cp_air * (e_sat_water_c / (pressure[0, i_c] - e_sat_water_c) - e_inf / (p_0 - e_inf))
	m_evap_c = m_evap_per_c * area[0, i_c]
	m_imp_c = m_imp[n_c, i_c]
	m_in_c = 0
	if (m_in[0, i_c - 1] + m_imp_c <= m_evap_c):
		m_in_c = 0
		m_evap_c = m_in_c + m_imp_c
		m_evap_per_c = m_evap_c / area[n_c, i_c]
	elif (T_wall_c + T0_c >= 368.15):
		m_in_c = 0
		m_evap_c = m_in[n_c, i_c - 1] + m_imp_c
		m_evap_per_c = m_evap_c / area[n_c, i_c]
	else:
		m_in_c = (m_in[0, i_c - 1] + m_imp_c - m_evap_c)
	m_out[n_c, i_c] = m_in_c
	Q_imp_per_c = m_imp_c * (u_inf * u_inf / 2 + Cp_water * T_inf)
	Q_imp_c = Q_imp_per_c * area[n_c, i_c]
	Q_evap_per_c = m_evap_per_c * (Le_c + Cp_water * T_s_c)
	Q_evap_c = Q_evap_per_c * area[n_c, i_c]
	Q_conv_per_c = h_air[n_c, i_c] * (T_s_c - 273.05)
	Q_conv_c = Q_conv_per_c * area[n_c, i_c]
	Q_anti_c = q_anti[n_c, i_c] * area[n_c, i_c]
	Q_in_c = 0

	T_ref_2 = (Q_in_c + Q_imp_c + Q_anti_c - Q_evap_c) / h_air[n_c, i_c] / area[n_c, i_c] + 247.8

	return T_ref_2










class MyProblem(ea.Problem):
	def __init__(self):
		name='MyProblem'  #初始化函数名称
		M=1  #初始化目标维数
		maxormins=[1]  #初始化最大最小值，1：最小化该目标；-1：最大化该目标
		Dim=4 #初始化决策变量维数
		varTypes=[0]*Dim #初始化决策变量类型，元素为0表示对应的变量是连续的；1表示是离散的，连续型变量
		lb=[10000 for x in range(0, 4)] #决策变量下界
		ub=[45000 for x in range(0, 4)]  #决策变量上界
		lbin=[0 for x in range(0,4)]  ## 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
		ubin=[1 for x in range(0, 4)]  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
		ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)













	def aimFunc(self, pop):

		Vars=pop.Phen#种群表现型，决策变量矩阵，每一列对应一个决策变量，每一行对应一个个体。
		#Vars=np.ones((NIND, 4))
		x = np.zeros((NIND, 4))#4区域加热， NIND种群规模
		for i in range(0, 4):
			x[:,[i]] = Vars[:, [i]]#每一行是一条染色体，也就是一个个体，一个具体问题的解，多列
		for i in range(150, 190):
			if ((i>=150) & (i<160)):
				q_anti[:,[i]]=x[:,[0]]
			if ((i>=160) & (i<170)):
				q_anti[:,[i]] = x[:,[1]]
			if ((i >= 170) & (i < 180)):
				q_anti[:,[i]]=x[:,[2]]
			if ((i >= 180) & (i <=190)):
				q_anti[:,[i]] = x[:,[3]]
		# for i in range(150, 154):
		# 	if ((i==150) ):
		# 		q_anti[:,[i]]=x[:,[0]]
		# 	if ((i==151) ):
		# 		q_anti[:,[i]] = x[:,[1]]
		# 	if ((i==152)):
		# 		q_anti[:,[i]]=x[:,[2]]
		# 	if ((i == 153)):
		# 		q_anti[:,[i]] = x[:,[3]]
		#for i in range(150,309):  #先试试4个网格
		for i in range(150, 190):
			if i==150:
				T_s[:,[i]] = 273.15
				T_wall[:,[i]] = T_s[:,[i]] - T0[:, [0]]
				for j in range(0, 10000):
					T_ref=0.5*(T_s[0, [i]]+T_wall[0, [i]]+T0[:, [0]]) #参考温度，此时参考温度是长为NIND的列向量
					for j_1 in range(0, NIND):
						if(T_ref[j_1, 0]>=T0[j_1,0]):#一个种群中的每一个个体都需要进行比较
							#rou_water[j_1,0]=(999.8396+18.224944*(T_ref[j_1, 0]-T0[j_1,0])-7.922210e-3*(T_ref[j_1, 0]-T0[j_1,0])**2.-55.44846e-6*(T_ref[j_1, 0]-T0[j_1,0])**3+149.7562e-9*(T_ref[j_1, 0]-T0[j_1,0])**4-393.2952e-12*(T_ref[j_1, 0]-T0[j_1,0])**5)/(1+(18.159725e-3)*(T_ref[j_1, 0]-T0[j_1,0]))
							rou_water[j_1,0]=998
							Cp_water[j_1,0]=4185.8518*(0.9979+3.1e-6*(T_ref[j_1, 0]-T0[j_1,0]-35)**2+3.8e-9*(T_ref[j_1, 0]-T0[j_1,0]-35)**4)
							miu_water[j_1,0]=1e-3*(1.76*m.exp(-3.5254e-2*(T_ref[j_1, 0]-T0[j_1,0])+4.7163e-4*(T_ref[j_1, 0]-T0[j_1,0])**2.-6.0667e-6*(T_ref[j_1, 0]-T0[j_1,0])**3))
						else:
							#rou_water[j_1,0]=1000 * (0.99986 + 6.69e-5 * (T_ref[j_1, 0]-T0[j_1,0]) - 8.486e-6 * (T_ref[j_1, 0]-T0[j_1,0]) ** 2.	+ 1.518e-7 * (T_ref[j_1, 0]-T0[j_1,0]) ** 3 - 6.9484e-9 * (T_ref[j_1, 0]-T0[j_1,0]) ** 4 - 3.6449e-10 ** (T_ref[j_1, 0]-T0[j_1,0]) ** 5.- 7.497e-12 * (T_ref[j_1, 0]-T0[j_1,0]) ** 6)
							rou_water[j_1, 0]=998
							Cp_water[j_1,0]=4185.8518 * (1.000938 - 2.7052e-3 * (T_ref[j_1, 0]-T0[j_1,0]) - 2.3235e-5 * (T_ref[j_1, 0]-T0[j_1,0]) ** 2.	+ 4.3778e-6 * (T_ref[j_1, 0]-T0[j_1,0]) ** 3 + 2.7136e-7 * (T_ref[j_1, 0]-T0[j_1,0]) ** 4)
							miu_water[j_1,0]=1e-3 * (1.76 * m.exp(-5.5721e-2 * (T_ref[j_1, 0]-T0[j_1,0]) - 1.3943e-3 * (T_ref[j_1, 0]-T0[j_1,0]) ** 2. - 4.3015e-5 * (T_ref[j_1, 0]-T0[j_1,0]) ** 3))
					for j_1_1 in range(0, NIND):
						e_sat_water[j_1_1,i] = 610.70 * m.exp(17.15 * (T_s[j_1_1,i] - T0[j_1_1,0]) / (T_s[j_1_1,i] - 38.25))
					Le[:,[i]]= 4185.8518 * (597.3 - 0.561 * (T_s[:,[i]] - T0[:,[0]]))
					m_evap_per[:,[i]] = 0.622 * h_air[:,[i]]/ Cp_air * (e_sat_water[:,[i]] / (pressure[:,[i]] - e_sat_water[:,[i]]) - e_inf / (p_0 - e_inf))
					m_evap[:,[i]]=m_evap_per[:,[i]]*area[:,[i]] 
					m_in[:,[i]]=0
					m_imp_per[:,[i]]=LWC * beta[:,[i]] * u_inf
					m_imp[:,[i]]=m_imp_per[:,[i]]*area[:,[i]]
					for j_2 in range(0, NIND):
						if(m_in[j_2, i]+m_imp[j_2, i]<=m_evap[j_2, i]):
							m_stat_p[j_2, i]=0
							m_evap[j_2, i]=m_in[j_2, i]+m_imp[j_2, i]
							m_evap_per[j_2, i]=m_evap[j_2, i]/area[j_2, i]
						elif (T_wall[j_2, i]+T0[j_2, 0]>=368.15):
							m_stat_p[j_2, i]=0
							m_evap[j_2, i]=m_in[j_2, i]+m_imp[j_2, i]
							m_evap_per[j_2, i] = m_evap[j_2, i] / area[j_2, i]
						else:
							m_stat_p[j_2, i]=m_imp[j_2, i]-m_evap[j_2, i]
						if m_stat_p[j_2, i]==0:
							delt_w[j_2, i]=0
					Q_in[:,[i]]=0
					Q_imp_per[:,[i]]=m_imp_per[:,[i]]*(u_inf*u_inf/2 + Cp_water[:,[0]]*T_inf)
					Q_imp[:,[i]]=Q_imp_per[:,[i]]*area[:,[i]]
					Q_evap_per[:,[i]]=m_evap_per[:,[i]]*(Le[:,[i]] +Cp_water[:,[0]]*T_s[:,[i]])
					Q_evap[:,[i]]=Q_evap_per[:,[i]]*area[:,[i]]
					Q_anti[:,[i]]=q_anti[:,[i]]*area[:,[i]]
					
					Q_conv_per[:,[i]]=h_air[:,[i]]*(T_s[:,[i]]-273.05)
					Q_conv[:,[i]]=Q_conv_per[:,[i]]*area[:,[i]]
					Q_sta_p[:,[i]]=(Q_imp[:,[i]]+Q_anti[:,[i]])-(Q_evap[:,[i]]+Q_conv[:,[i]])
					T_ref_1[:,[0]] = Q_sta_p[:,[i]] / m_stat_p[:,[i]] / Cp_water[:,[0]]
					for j_3 in range(0, NIND):
						while (T_s[j_3, i]-T_ref_1[j_3, 0]>0.00001):
							 T_s[j_3, i]=T_s[j_3, i]+0.01*(T_ref_1[j_3, 0]-T_s[j_3, i])
					err[:,[0]]=(abs(T_s[:,[i]]-T_ref_1[:, [0]]))
					if (((err[:, [0]])<=0.00001).all()):
						break
					else:
						T_s[:,[i]]=T_s[:,[i]]+0.01*(T_ref_1[:, [0]]-T_s[:,[i]])
				T_wall[:,[i]]=T_s[:,[i]]-T0[:,[0]]
				# print(T_wall[-1, 150])
			else:
				T_s[:,[i]] = T_s[:, [i-1]]
				# for f in range(0, NIND):
				# 	if (T_s[f,i]>400):
				# 		T_s[f, i]=310
				# 	elif(T_s[f, i]<200):
				# 		T_s[f, i]=310
				# print(T_s[:, [i]])

				T_wall[:,[i]]=T_wall[:, [i-1]]
				# print("T_wall:%f"%(T_wall[-1, i]))
				# print("计算到第 %d 个网格" % (i))
				for k in range(0, 100):
					T_ref=0.5*(T_s[0, [i]]+T_wall[0, [i]]+T0[:,[0]])
					# print(T_ref)
					m_imp_per[:,[i]]=LWC * beta[:,[i]] * u_inf
					m_imp[:,[i]]=m_imp_per[:,[i]]*area[:,[i]]

					for k_1 in range(0, NIND):
						if(T_ref[k_1, 0]>=T0[k_1,0]):#一个种群中的每一个个体都需要进行比较
							#rou_water[k_1,0]=(999.8396+18.224944*(T_ref[k_1, 0]-T0[k_1,0])-7.922210e-3*(T_ref[k_1, 0]-T0[k_1,0])**2.-55.44846e-6*(T_ref[k_1, 0]-T0[k_1, 0])**3+149.7562e-9*(T_ref[k_1, 0]-T0[k_1, 0])**4-393.2952e-12*(T_ref[k_1, 0]-T0[k_1, 0])**5)/(1+(18.159725e-3)*(T_ref[k_1, 0]-T0[k_1, 0]))
							rou_water[k_1,0]=998
							Cp_water[k_1,0]=4185.8518*(0.9979+3.1e-6*(T_ref[k_1,0]-T0[k_1,0]-35)**2+3.8e-9*(T_ref[k_1, 0]-T0[k_1,0]-35)**4)
							miu_water[k_1,0]=1e-3*(1.76*m.exp(-3.5254e-2*(T_ref[k_1, 0]-T0[k_1,0])+4.7163e-4*(T_ref[k_1, 0]-T0[k_1,0])**2.-6.0667e-6*(T_ref[k_1, 0]-T0[k_1,0])**3))
						else:
							#rou_water[k_1,0]=1000 * (0.99986 + 6.69e-5 * (T_ref[k_1, 0]-T0[k_1,0]) - 8.486e-6 * (T_ref[k_1, 0]-T0[k_1,0]) ** 2.	+ 1.518e-7 * (T_ref[k_1, 0]-T0[k_1,0]) ** 3 - 6.9484e-9 * (T_ref[k_1, 0]-T0[k_1,0]) ** 4 - 3.6449e-10 ** (T_ref[k_1, 0]-T0[k_1,0]) ** 5.- 7.497e-12 * (T_ref[k_1, 0]-T0[k_1,0]) ** 6)
							rou_water[k_1,0]=998
							Cp_water[k_1,0]=4185.8518 * (1.000938 - 2.7052e-3 * (T_ref[k_1, 0]-T0[k_1,0]) - 2.3235e-5 * (T_ref[k_1, 0]-T0[k_1,0]) ** 2.	+ 4.3778e-6 * (T_ref[k_1, 0]-T0[k_1,0]) ** 3 + 2.7136e-7 * (T_ref[k_1, 0]-T0[k_1,0]) ** 4)
							miu_water[k_1,0]=1e-3 * (1.76 * m.exp(-5.5721e-2 * (T_ref[k_1, 0]-T0[k_1,0]) - 1.3943e-3 * (T_ref[k_1, 0]-T0[k_1,0]) ** 2. - 4.3015e-5 * (T_ref[k_1, 0]-T0[k_1,0]) ** 3))

						e_sat_water[k_1,i] = 610.70 * m.exp(17.15 * (T_s[k_1,i] - T0[k_1, 0]) / (T_s[k_1,i] - 38.25))
						# print("NIND: %d"%(k_1))
						# print("Q_anti: %f" % (Q_anti[k_1, i-1]))
						# print("Q_in: %f" % (Q_in[k_1, i-1]))
						# print("e_inf: %f" %(e_inf))
						# print("p_0: %f" %(p_0))
						# print("e_sat_water= %f" % (e_sat_water[k_1, i-1]))
						#
						# print("m_evap: %f" % (m_evap[k_1, i - 1]))
						# print("m_evap-1: %f" % (m_evap[k_1, i - 2]))
						# print("Q_evap: %f" % (Q_evap[k_1, i-1]))
						# print("Q_evap-1: %f" % (Q_evap[k_1, i - 2]))
						# print("Q_imp: %f" % (Q_imp[k_1, i - 1]))
						#
						# print("T_ref2: %f" % (T_ref_2[k_1, 0]))
						# print("T_s=%f" %(T_s[k_1, i]))
						# print("P_s=%f" %(q_anti[k_1, i-1]) )
						#
						# print("计算到第 %d 个网格" % (i))

					Le[:,[i]]= 4185.8518 * (597.3 - 0.561 * (T_s[:,[i]] - T0[:,[0]]))

					m_evap_per[:,[i]] = 0.622 * h_air[:,[i]]/ Cp_air * ((e_sat_water[:,[i]] / (pressure[:,[i]] - e_sat_water[:,[i]])) - e_inf / (p_0 - e_inf))

					m_evap[:,[i]]=m_evap_per[:,[i]]*area[:,[i]]
					for k_2 in range(0, NIND):
						if(m_in[k_2, i]+m_imp[k_2, i]<=m_evap[k_2, i]):
							m_in[k_2, i]=0
							m_evap[k_2, i]=m_in[k_2, i]+m_imp[k_2, i]
							m_evap_per[k_2, i]=m_evap[k_2, i]/area[k_2, i]
						elif (T_wall[k_2, i]+T0[k_2, 0]>=368.15):
							m_in[k_2, i]=0
							m_evap[k_2, i]=m_in[k_2, i]+m_imp[k_2, i]
							m_evap_per[k_2, i] = m_evap[k_2, i] / area[k_2, i]
						else:
							m_in[k_2, i]=(m_in[k_2, i-1]+m_imp[k_2, i]-m_evap[k_2, i])

					m_out[:,[i-1]]=m_in[:,[i]]
					# print("m_out: %f" %(m_out[-1, i]))
					Q_imp_per[:,[i]]=m_imp_per[:,[i]]*(u_inf*u_inf/2 + Cp_water[:,[0]]*T_inf)
					Q_imp[:,[i]]=Q_imp_per[:,[i]]*area[:,[i]]
					Q_evap_per[:,[i]]=m_evap_per[:,[i]]*(Le[:,[i]] +Cp_water[:,[0]]*T_s[:,[i]])
					Q_evap[:,[i]]=Q_evap_per[:,[i]]*area[:,[i]]
					
					Q_conv_per[:,[i]]=h_air[:,[i]]*(T_s[:,[i]]-273.05)
					Q_conv[:,[i]]=Q_conv_per[:,[i]]*area[:,[i]]
					Q_anti[:,[i]]=q_anti[:,[i]]*area[:,[i]]
					
					T_ref_2[:, [0]]=(Q_in[:, [i-1]]+Q_imp[:, [i]]+Q_anti[:, [i]] -Q_evap[:, [i]])/h_air[:, [i]]/area[:, [i]] + 247.8

					for k_3 in range(0, NIND):
						while (abs(T_s[k_3, i] - T_ref_2[k_3, 0]) > 0.00001):
							T_ref_2[k_3, 0]=cell_iteration(i, k_3, T_wall[k_3, i], T_s[k_3, i])  #这个的参数必须仔细斟酌。
							# print("T_ref_2: %f, NIND_num: %d" %(T_ref_2[k_3, 0], k_3) )
							# print("T_s: %f, NIND_num: %d, 第 %d 个网格, 第 %d 次迭代" % (T_s[k_3, i], k_3, i ,k))

							T_s[k_3 , i]=T_s[k_3 , i]+0.01 * (T_ref_2[k_3, 0]-T_s[k_3 , i])

						T_wall[k_3, i]=T_s[k_3 , i]-273.15


				
				#------begin这部分个体不收敛的迭代过程是错误的。这个过程只是单独+误差，并没有将误差参与迭代过程。
				# for k_3 in range(0, NIND):
				# 	while (abs(T_s[k_3, i] - T_ref_2[k_3, 0] )> 0.00001):
				# 		T_s[k_3, i] = T_s[k_3, i] + 0.01 * (T_ref_2[k_3, 0] - T_s[k_3, i])
				# err[:, [0]] = (abs(T_s[:, [i]] - T_ref_2[:, [0]]))
				# if (((err[:, [0]]) <= 0.0001).all()):
				# 	break
				# else:
				# 	T_s[:, [i]] = T_s[:, [i]] + 0.01 * (T_ref_1[:, [0]] - T_s[:, [i]])
				# ----------end----------------------------------------------------

				#---------begin先找出没有收敛的个体，然后+误差，最后整个种群迭代。
				# for k_3 in range(0, NIND):
				# 	if(abs(T_s[k_3, i] - T_ref_2[k_3, 0] )> 0.00001):
				# 		T_s[k_3, i] = T_s[k_3, i] + 0.01 * (T_ref_2[k_3, 0] - T_s[k_3, i])
				#
				# err[:, [0]] = (abs(T_s[:, [i]] - T_ref_2[:, [0]]))
				# if (((err[:, [0]]) <= 0.0001).all()):
				# 	break
				# else:
				# 	T_s[:, [i]] = T_s[:, [i]] + 0.01 * (T_ref_1[:, [0]] - T_s[:, [i]])
				#----------end----------------------------------------------------








					#---------begin修改后的个体迭代过程------------
					# err[:, [0]] = (abs(T_s[:, [i]] - T_ref_2[:, [0]]))
					# if (((err[:, [0]]) <= 0.0001).all()):
					# 	# print("T_s(k_3)=%f" % (T_s[-1]))
					# 	break
					# else:
					# 	for k_4 in range(0, NIND):
					# 		if(abs(T_s[k_4, i] - T_ref_2[k_4, 0])>0.0001):
					# 			T_s[k_4,i]=T_s[k_4,i]+0.01*(T_ref_2[k_4, 0]-T_s[k_4,i])
					#----------------end修改后的个体迭代过程----------------






				T_wall[:, [i]] = T_s[:, [i]] - T0[:, [0]]
			print('T_wall_151',T_wall[-1, 151])
			print('T_wall_189', T_wall[-1, 189])
			print('m_in_189:', m_in[-1, 189] )
			print('A区功率:%f' %(q_anti[-1, 151]))
			print('D区功率:%f' % (q_anti[-1, 189]))
		

		aim2=np.zeros((NIND, 1))
		for i in range(150, 190):
			aim2[:,[0]]+=q_anti[:,[i]]*area[:, [i]]
		pop.ObjV=aim2[:,[0]]
		
		pop.CV=np.hstack([
							# 25-T_wall[:, [189]],
							23-T_wall[:, [151]],
							0.0003-m_in[:, [189]]

						])


			

