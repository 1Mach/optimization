import numpy as np
import geatpy as ea
import math as m

n=309
NIND=20
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

x = np.zeros((NIND, 4))  # 4区域加热， NIND种群规模
x[:, [0]]=20.15e3
x[:,[1]]=21.7e3
for i in range(150, 152):
	if (i == 150) :
		q_anti[:, [i]] = x[:, [0]]
	if (i == 151) :
		q_anti[:, [i]] = x[:, [1]]

for i in range(150,152):
			if i==150:
				T_s[:,[i]] = 273.15
				T_wall[:,[i]] = T_s[:,[i]] - T0[:, [0]]
				for j in range(0, 100):
					T_ref=0.5*(T_s[0, [i]]+T_wall[0, [i]]+T0[:, [0]]) #参考温度，此时参考温度是长为NIND的列向量
					for j_1 in range(0, NIND):
						if(T_ref[j_1, 0]>=T0[j_1,0]):#一个种群中的每一个个体都需要进行比较
							rou_water[j_1,0]=(999.8396+18.224944*(T_ref[j_1, 0]-T0[j_1,0])-7.922210e-3*(T_ref[j_1, 0]-T0[j_1,0])**2.-55.44846e-6*(T_ref[j_1, 0]-T0[j_1,0])**3+149.7562e-9*(T_ref[j_1, 0]-T0[j_1,0])**4-393.2952e-12*(T_ref[j_1, 0]-T0[j_1,0])**5)/(1+(18.159725e-3)*(T_ref[j_1, 0]-T0[j_1,0]))
							Cp_water[j_1,0]=4185.8518*(0.9979+3.1e-6*(T_ref[j_1, 0]-T0[j_1,0]-35)**2+3.8e-9*(T_ref[j_1, 0]-T0[j_1,0]-35)**4)
							miu_water[j_1,0]=1e-3*(1.76*m.exp(-3.5254e-2*(T_ref[j_1, 0]-T0[j_1,0])+4.7163e-4*(T_ref[j_1, 0]-T0[j_1,0])**2.-6.0667e-6*(T_ref[j_1, 0]-T0[j_1,0])**3))
						else:
							rou_water[j_1,0]=1000 * (0.99986 + 6.69e-5 * (T_ref[j_1, 0]-T0[j_1,0]) - 8.486e-6 * (T_ref[j_1, 0]-T0[j_1,0]) ** 2.	+ 1.518e-7 * (T_ref[j_1, 0]-T0[j_1,0]) ** 3 - 6.9484e-9 * (T_ref[j_1, 0]-T0[j_1,0]) ** 4 - 3.6449e-10 ** (T_ref[j_1, 0]-T0[j_1,0]) ** 5.- 7.497e-12 * (T_ref[j_1, 0]-T0[j_1,0]) ** 6)
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
				print(T_wall[:,[i]])
			else:
				T_s[:,[i]] = T_s[:, [i-1]]
				T_wall[:,[i]]=T_wall[:, [i-1]]
				for k in range(0, 100):
					T_ref=0.5*(T_s[0, [i]]+T_wall[0, [i]]+T0[:,[0]])
					m_imp_per[:,[i]]=LWC * beta[:,[i]] * u_inf
					m_imp[:,[i]]=m_imp_per[:,[i]]*area[:,[i]]
					for k_1 in range(0, NIND):
						if(T_ref[k_1, 0]>=T0[k_1,0]):#一个种群中的每一个个体都需要进行比较
							rou_water[k_1,0]=(999.8396+18.224944*(T_ref[k_1, 0]-T0[k_1,0])-7.922210e-3*(T_ref[k_1, 0]-T0[k_1,0])**2.-55.44846e-6*(T_ref[k_1, 0]-T0[k_1, 0])**3+149.7562e-9*(T_ref[k_1, 0]-T0[k_1, 0])**4-393.2952e-12*(T_ref[k_1, 0]-T0[k_1, 0])**5)/(1+(18.159725e-3)*(T_ref[k_1, 0]-T0[k_1, 0]))
							Cp_water[k_1,0]=4185.8518*(0.9979+3.1e-6*(T_ref[k_1,0]-T0[k_1,0]-35)**2+3.8e-9*(T_ref[k_1, 0]-T0[k_1,0]-35)**4)
							miu_water[k_1,0]=1e-3*(1.76*m.exp(-3.5254e-2*(T_ref[k_1, 0]-T0[k_1,0])+4.7163e-4*(T_ref[k_1, 0]-T0[k_1,0])**2.-6.0667e-6*(T_ref[k_1, 0]-T0[k_1,0])**3))
						else:
							rou_water[k_1,0]=1000 * (0.99986 + 6.69e-5 * (T_ref[k_1, 0]-T0[k_1,0]) - 8.486e-6 * (T_ref[k_1, 0]-T0[k_1,0]) ** 2.	+ 1.518e-7 * (T_ref[k_1, 0]-T0[k_1,0]) ** 3 - 6.9484e-9 * (T_ref[k_1, 0]-T0[k_1,0]) ** 4 - 3.6449e-10 ** (T_ref[k_1, 0]-T0[k_1,0]) ** 5.- 7.497e-12 * (T_ref[k_1, 0]-T0[k_1,0]) ** 6)
							Cp_water[k_1,0]=4185.8518 * (1.000938 - 2.7052e-3 * (T_ref[k_1, 0]-T0[k_1,0]) - 2.3235e-5 * (T_ref[k_1, 0]-T0[k_1,0]) ** 2.	+ 4.3778e-6 * (T_ref[k_1, 0]-T0[k_1,0]) ** 3 + 2.7136e-7 * (T_ref[k_1, 0]-T0[k_1,0]) ** 4)
							miu_water[k_1,0]=1e-3 * (1.76 * m.exp(-5.5721e-2 * (T_ref[k_1, 0]-T0[k_1,0]) - 1.3943e-3 * (T_ref[k_1, 0]-T0[k_1,0]) ** 2. - 4.3015e-5 * (T_ref[k_1, 0]-T0[k_1,0]) ** 3))
					for k_1_1 in range(0, NIND):
						e_sat_water[k_1_1,i] = 610.70 * m.exp(17.15 * (T_s[k_1_1,i] - T0[k_1_1, 0]) / (T_s[k_1_1,i] - 38.25))
					Le[:,[i]]= 4185.8518 * (597.3 - 0.561 * (T_s[:,[i]] - T0[:,[0]]))
					m_evap_per[:,[i]] = 0.622 * h_air[:,[i]]/ Cp_air * (e_sat_water[:,[i]] / (pressure[:,[i]] - e_sat_water[:,[i]]) - e_inf / (p_0 - e_inf))
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
					m_out[:,[i-1]]=m_in[:,[i-1]]

					Q_imp_per[:,[i]]=m_imp_per[:,[i]]*(u_inf*u_inf/2 + Cp_water[:,[0]]*T_inf)
					Q_imp[:,[i]]=Q_imp_per[:,[i]]*area[:,[i]]
					Q_evap_per[:,[i]]=m_evap_per[:,[i]]*(Le[:,[i]] +Cp_water[:,[0]]*T_s[:,[i]])
					Q_evap[:,[i]]=Q_evap_per[:,[i]]*area[:,[i]]

					Q_conv_per[:,[i]]=h_air[:,[i]]*(T_s[:,[i]]-273.05)
					Q_conv[:,[i]]=Q_conv_per[:,[i]]*area[:,[i]]
					Q_anti[:,[i]]=q_anti[:,[i]]*area[:,[i]]

					T_ref_2[:, [0]]=(Q_in[:, [i-1]]+Q_imp[:, [i]]+Q_anti[:, [i]] -Q_evap[:, [i]])/h_air[:, [i]]/area[:, [i]] + 247.8
					for k_3 in range(0, NIND):
						while (T_s[k_3, i]-T_ref_2[k_3, 0]>0.00001):
							 T_s[k_3, i]=T_s[k_3, i]+0.01*(T_ref_2[k_3, 0]-T_s[k_3, i])
					err_2=(abs(T_s[:,[i]]-T_ref_2[:, [0]]))
					if (((err_2[:, [0]])<=0.00001).all()):
						break
					else:
						T_s[:,[i]]=T_s[:,[i]]+0.01*(T_ref_2[:, [0]]-T_s[:,[i]])
aim2=0
for i in range(150, 309):
	aim2 += q_anti[:, [i]] * area[:, [i]]
