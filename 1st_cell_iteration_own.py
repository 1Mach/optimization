#利用自己的计算结果进行第一个网格的迭代过程

import numpy as np
import geatpy as ea
import math as m

#-----------------空气参数----------------------
R_a           = 287
Cp_air        = 1013



#-----------------------实验参数---------------------------
LWC           = 0.55*1e-3
MVD           = 20*1e-6
p_tot         = 70000
u_inf         = 89.4
T_tot         = 251.55
T0            = 273.15
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

#液膜表面饱和水蒸气压
e_inf=610.70*m.exp(17.15*(T_inf-T0)/(T_inf-38.25))



n=309



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
#area=np.repeat(area, NIND, axis=0)#遗传算法矩阵复制
# print(area)
betalist.append(read_original_data('./original_data/beta'))
beta=np.asarray(betalist)
#beta=np.repeat(beta, NIND, axis=0)
h_airlist.append(read_original_data('./original_data/heat_transfer_coefficient'))
h_air=np.asarray(h_airlist)
#h_air=np.repeat(h_air, NIND, axis=0)
pressurelist.append(read_original_data('./original_data/pressure'))
pressure=np.asarray(pressurelist)
#pressure=np.repeat(pressure, NIND, axis=0)


Twall_steadylist.append(read_original_data('./original_data/temperature'))
# print(Twall_steadylist)
Twall_steady=np.asarray(Twall_steadylist)#计算的稳态的温度
# print(Twall_steady[0,3])
T_wall=np.zeros((1, n))

T_s=np.zeros((1, n))

Le=np.zeros((1, n))
e_sat_water=np.zeros((1, n))
m_evap_per=np.zeros((1, n))
m_evap=np.zeros((1, n))
m_in=np.zeros((1, n))
m_out=np.zeros((1, n))
m_stat_p=np.zeros((1, n))
m_imp=np.zeros((1, n))
m_imp_per=np.zeros((1, n))
delt_w=np.zeros((1, n))
Q_in=np.zeros((1, n))
Q_out=np.zeros((1, n))
Q_imp=np.zeros((1, n))
Q_imp_per=np.zeros((1, n))
Q_evap=np.zeros((1, n))
Q_evap_per=np.zeros((1, n))
Q_anti=np.zeros((1, n))
Q_conv=np.zeros((1, n))
Q_conv_per=np.zeros((1, n))
Q_sta_p=np.zeros((1, n))
q_anti=np.zeros((1, n))

i=150
T_s[0, i]=273.15
T_wall[0, i]=T_s[0, i]-T0
for i in range(150, 309):
	if ((i>=150) & (i<=180)):
		q_anti[0, i]=20.15*1e+3
	if ((i>180) & (i<=210)):
		q_anti[0, i] = 26.35 * 1e+3
	if ((i > 210) & (i <= 260)):
		q_anti[0, i]=18.6* 1e+3
	if ((i > 260) & (i < 309)):
		q_anti[0, i] = 18.6 * 1e+3



for i in range(150, 309):
	if i==150:
		T_s[0, i] = 273.15
		T_wall[0, i] = T_s[0, i] - T0
		q_anti[0, i] = 20.15 * 1e+3
		for j in range(0, 100000):
			T_ref=0.5*(T_s[0, i]+T_wall[0, i]+T0) #参考温度
			if(T_ref>=T0):
				rou_water=(999.8396+18.224944*(T_ref-T0)-7.922210e-3*(T_ref-T0)**2.-55.44846e-6*(T_ref-T0)**3+149.7562e-9*(T_ref-T0)**4-393.2952e-12*(T_ref-T0)**5)/(1+(18.159725e-3)*(T_ref-T0))
				Cp_water=4185.8518*(0.9979+3.1e-6*(T_ref-T0-35)**2+3.8e-9*(T_ref-T0-35)**4)
				miu_water=1e-3*(1.76*m.exp(-3.5254e-2*(T_ref-T0)+4.7163e-4*(T_ref-T0)**2.-6.0667e-6*(T_ref-T0)**3))
			else:
				rou_water = 1000 * (0.99986 + 6.69e-5 * (T_ref - T0) - 8.486e-6 * (T_ref - T0) ** 2.	+ 1.518e-7 * (T_ref - T0) ** 3 - 6.9484e-9 * (T_ref - T0) ** 4 - 3.6449e-10 ** (T_ref - T0) ** 5.- 7.497e-12 * (T_ref - T0) ** 6)
				Cp_water = 4185.8518 * (1.000938 - 2.7052e-3 * (T_ref - T0) - 2.3235e-5 * (T_ref - T0) ** 2.	+ 4.3778e-6 * (T_ref - T0) ** 3 + 2.7136e-7 * (T_ref - T0) ** 4)
				miu_water = 1e-3 * (1.76 * m.exp(-5.5721e-2 * (T_ref - T0) - 1.3943e-3 * (T_ref - T0) ** 2. - 4.3015e-5 * (T_ref - T0) ** 3))
			#饱和蒸汽压
			e_sat_water[0, i ] = 610.70 * m.exp(17.15 * (T_s[0, i] - T0) / (T_s[0, i] - 38.25))

			#蒸发热
			Le[0, i]= 4185.8518 * (597.3 - 0.561 * (T_s[0, i] - T0))
			m_evap_per[0, i] = 0.622 * h_air[0, i]/ Cp_air * (e_sat_water[0, i] / (pressure[0, i] - e_sat_water[0, i]) - e_inf / (p_0 - e_inf))
			m_evap[0, i]=m_evap_per[0, i]*area[0, i] #这个需要替换成我自己的结果
			m_in[0, i]=0
			m_imp_per[0, i]=LWC * beta[0, i] * u_inf
			m_imp[0, i]=m_imp_per[0, i]*area[0, i]
			if (m_in[0, i]+m_imp[0, i]<=m_evap[0, i]):
				m_stat_p[0, i]=0
				m_evap[0, i]=m_in[0, i]+m_imp[0, i]
				m_evap_per[0, i]=m_evap[0, i]/area[0, i]
			elif(T_wall[0, i]+T0>=368.15):
				m_stat_p[0, i]=0
				m_evap[0, i]=m_in[0, i]+m_imp[0, i]
				m_evap_per[0, i] = m_evap[0, i] / area[0, i]
			else:
				m_stat_p[0, i]=m_imp[0, i]-m_evap[0, i]

			if m_stat_p[0, i]==0:
				delt_w[0, i]=0
			Q_in[0, i]=0
			Q_imp_per[0, i]=m_imp_per[0, i]*(u_inf*u_inf/2 + Cp_water*T_inf)
			Q_imp[0, i]=Q_imp_per[0, i]*area[0, i]
			Q_evap_per[0, i]=m_evap_per[0, i]*(Le[0, i] +Cp_water*T_s[0, i])
			Q_evap[0, i]=Q_evap_per[0, i]*area[0, i]
			Q_anti[0, i]=q_anti[0, i]*area[0, i]

			Q_conv_per[0, i]=h_air[0, i]*(T_s[0, i]-273.05)
			Q_conv[0, i]=Q_conv_per[0, i]*area[0, i]
			Q_sta_p[0, i]=(Q_imp[0, i]+Q_anti[0, i])-(Q_evap[0, i]+Q_conv[0, i])

		#	T_ref_1=Q_sta_p[0, i]/m_stat_p[0, i]/Cp_water
			T_ref_1 = Q_sta_p[0, i] / m_stat_p[0, i] / Cp_water
			if ((abs(T_s[0, i]-T_ref_1))<=0.00001):
				break
			else:
				T_s[0, i]=T_s[0, i]+0.01*(T_ref_1-T_s[0, i])
		T_wall[0, i]=T_s[0, i]-T0
	else:
		T_s[0, i]=T_s[0, i-1]
		T_wall[0, i]=T_wall[0, i-1]
		for j_1 in range(0, 100000):
			T_ref=0.5*(T_s[0, i]+T_wall[0, i]+T0)
			m_imp_per[0, i] = LWC * beta[0, i] * u_inf
			m_imp[0, i] = m_imp_per[0, i] * area[0, i]
			if T_ref>=T0:
				rou_water=(999.8396+18.224944*(T_ref-T0)-7.922210e-3*(T_ref-T0)**2 -55.44846e-6*(T_ref-T0)**3+149.7562e-9*(T_ref-T0)**4-393.2952e-12*(T_ref-T0)**5)/(1+(18.159725e-3)*(T_ref-T0))
				Cp_water=4185.8518*(0.9979+3.1e-6*(T_ref-T0-35)**2+3.8e-9*(T_ref-T0-35)**4)
				miu_water=1e-3*(1.76*m.exp(-3.5254e-2*(T_ref-T0)+4.7163e-4*(T_ref-T0)**2.-6.0667e-6*(T_ref-T0)**3))
			else:
				rou_water = 1000 * (0.99986 + 6.69e-5 * (T_ref - T0) - 8.486e-6 * (T_ref - T0) ** 2.	+ 1.518e-7 * (T_ref - T0) ** 3 - 6.9484e-9 * (T_ref - T0) ** 4 - 3.6449e-10 ** (T_ref - T0) ** 5.- 7.497e-12 * (T_ref - T0) ** 6)
				Cp_water = 4185.8518 * (1.000938 - 2.7052e-3 * (T_ref - T0) - 2.3235e-5 * (T_ref - T0) ** 2.	+ 4.3778e-6 * (T_ref - T0) ** 3 + 2.7136e-7 * (T_ref - T0) ** 4)
				miu_water = 1e-3 * (1.76 * m.exp(-5.5721e-2 * (T_ref - T0) - 1.3943e-3 * (T_ref - T0) ** 2. - 4.3015e-5 * (T_ref - T0) ** 3))
			e_sat_water[0, i] = 610.70 * m.exp(17.15 * (T_s[0, i] - T0) / (T_s[0, i] - 38.25))
			Le[0, i] = 4185.8518 * (597.3 - 0.561 * (T_s[0, i] - T0))
			m_evap_per[0, i] = 0.622 * h_air[0, i]/ Cp_air * (e_sat_water[0, i] / (pressure[0, i] - e_sat_water[0, i]) - e_inf / (p_0 - e_inf))
			m_evap[0, i]=m_evap_per[0, i]*area[0, i]
			if (m_in[0, i-1]+m_imp[0, i]<=m_evap[0, i]):
				m_in[0, i]=0
				m_evap[0, i]=m_in[0, i]+m_imp[0, i]
				m_evap_per[0, i]=m_evap[0, i]/area[0, i]
			elif(T_wall[0, i]+T0>=368.15):
				m_in[0, i]=0
				m_evap[0, i]=m_in[0, i-1]+m_imp[0, i]
				m_evap_per[0, i] = m_evap[0, i] / area[0, i]
			else:
				m_in[0, i]=(m_in[0, i-1]+m_imp[0, i]-m_evap[0, i])
			m_out[0,i-1]=m_in[0, i]

			Q_imp_per[0, i]=m_imp_per[0, i]*(u_inf*u_inf/2 + Cp_water*T_inf)
			Q_imp[0, i]=Q_imp_per[0, i]*area[0, i]
			Q_evap_per[0, i]=m_evap_per[0, i]*(Le[0, i] +Cp_water*T_s[0, i])
			Q_evap[0, i]=Q_evap_per[0, i]*area[0, i]
			Q_conv_per[0, i]=h_air[0, i]*(T_s[0, i]-273.05)
			Q_conv[0, i]=Q_conv_per[0, i]*area[0, i]
			Q_anti[0, i]=q_anti[0, i]*area[0, i]

			# Q_in[0, i]=Q_in[0, i-1]+(Q_imp[0, i]+Q_anti[0, i])-(Q_evap[0, i]+Q_conv[0, i])
			# Q_out[0, i-1]=Q_in[0, i]

			T_ref_2=(Q_in[0, i-1]+Q_imp[0, i]+Q_anti[0, i] -Q_evap[0, i])/h_air[0, i]/area[0, i] + 247.8
			if (abs(T_s[0, i]-T_ref_2)<0.00001):
				break
			else:
				T_s[0, i]=T_s[0, i]+0.01 * (T_ref_2-T_s[0, i])
		T_wall[0, i]=T_s[0, i]-T0
print(e_sat_water)
print(T_s[0, 151])












