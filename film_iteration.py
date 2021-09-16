import numpy as np
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
# area=np.repeat(area, NIND, axis=0)#复制行

betalist.append(read_original_data('./original_data/beta'))
beta=np.asarray(betalist)
# beta=np.repeat(beta, NIND, axis=0)

h_airlist.append(read_original_data('./original_data/heat_transfer_coefficient'))
h_air=np.asarray(h_airlist)
# h_air=np.repeat(h_air, NIND, axis=0)


pressurelist.append(read_original_data('./original_data/pressure'))
pressure=np.asarray(pressurelist)
# pressure=np.repeat(pressure, NIND, axis=0)
T0=273.15
T_wall=np.zeros((1, n))

T_s=np.zeros((1, n))

# Le=np.zeros((1, n))
# e_sat_water=np.zeros((1, n))
# m_evap_per=np.zeros((1, n))
# m_evap=np.zeros((1, n))
m_in=np.zeros((1, n))
m_out=np.zeros((1, n))
# m_stat_p=np.zeros((1, n))
m_imp=np.zeros((1, n))
m_imp_per=np.zeros((1, n))
# delt_w=np.zeros((1, n))
# Q_in=np.zeros((1, n))
# Q_out=np.zeros((1, n))
# Q_imp=np.zeros((1, n))
# Q_imp_per=np.zeros((1, n))
# Q_evap=np.zeros((1, n))
# Q_evap_per=np.zeros((1, n))
# Q_anti=np.zeros((1, n))
# Q_conv=np.zeros((1, n))
# Q_conv_per=np.zeros((1, n))
# Q_sta_p=np.zeros((1, n))
# q_anti=np.zeros((1, n))


def cell_iteration(cell_num, T_s): #在一个网格内不断迭代收敛的函数。可以用在同一个网格的不同种群导致的T_s和T_ref不收敛时使用。
	i=cell_num

	T_ref = 0.5 * (T_s + T_wall_c + T0)

	if (T_ref >= T0):
		#rou_water=(999.8396+18.224944*(T_ref-T0)-7.922210e-3*(T_ref-T0)**2.-55.44846e-6*(T_ref-T0)**3+149.7562e-9*(T_ref-T0)**4-393.2952e-12*(T_ref-T0)**5)/(1+(18.159725e-3)*(T_ref-T0))
		Cp_water=4185.8518*(0.9979+3.1e-6*(T_ref-T0-35)**2+3.8e-9*(T_ref-T0-35)**4)
		# miu_water=1e-3*(1.76*m.exp(-3.5254e-2*(T_ref-T0)+4.7163e-4*(T_ref-T0)**2.-6.0667e-6*(T_ref-T0)**3))
	else:
		# rou_water = 1000 * (0.99986 + 6.69e-5 * (T_ref - T0) - 8.486e-6 * (T_ref - T0) ** 2.	+ 1.518e-7 * (T_ref - T0) ** 3 - 6.9484e-9 * (T_ref - T0) ** 4 - 3.6449e-10 ** (T_ref - T0) ** 5.- 7.497e-12 * (T_ref - T0) ** 6)
		Cp_water = 4185.8518 * (1.000938 - 2.7052e-3 * (T_ref - T0) - 2.3235e-5 * (T_ref - T0) ** 2.	+ 4.3778e-6 * (T_ref - T0) ** 3 + 2.7136e-7 * (T_ref - T0) ** 4)
		# miu_water = 1e-3 * (1.76 * m.exp(-5.5721e-2 * (T_ref - T0) - 1.3943e-3 * (T_ref - T0) ** 2. - 4.3015e-5 * (T_ref - T0) ** 3))
		#饱和蒸汽压
	e_sat_water = 610.70 * m.exp(17.15 * (T_s - T0) / (T_s - 38.25))
	Le = 4185.8518 * (597.3 - 0.561 * (T_s - T0))
	m_evap_per = 0.622 * h_air[0, i]/ Cp_air * (e_sat_water/ (pressure[0, i] - e_sat_water) - e_inf / (p_0 - e_inf))
	m_evap=m_evap_per*area[0, i]
	m_imp_c=m_imp[0, i]
	m_in_c=0
	if (m_in[0, i-1]+m_imp_c<=m_evap):
		m_in_c=0
		m_evap=m_in_c+m_imp_c
		m_evap_per=m_evap/area[0, i]
	elif(T_wall_c+T0>=368.15):
		m_in_c=0
		m_evap=m_in[0, i-1]+m_imp_c
		m_evap_per = m_evap[0, i] / area[0, i]
	else:
		m_in_c=(m_in[0, i-1]+m_imp[0, i]-m_evap)
	m_out[0,i-1]=m_in_c

	Q_imp_per=m_imp_c*(u_inf*u_inf/2 + Cp_water*T_inf)
	Q_imp=Q_imp_per*area[0, i]
	Q_evap_per=m_evap_per*(Le +Cp_water*T_s)
	Q_evap=Q_evap_per*area[0, i]
	Q_conv_per=h_air[0, i]*(T_s-273.05)
	Q_conv=Q_conv_per*area[0, i]
	q_anti=20150
	Q_anti=q_anti*area[0, i]
	Q_in_0=0


	T_ref_2 = (Q_in_0 + Q_imp+ Q_anti- Q_evap) / h_air[0, i] / area[0, i] + 247.8

	return T_ref_2


def cell_iteration_c(cell_num, NIND_num, T_wall_c, T_s_c):
	i_c = cell_num
	n_c=NIND_num

	T_ref_c = 0.5 * (T_s_c + T_wall_c + T0)
	if(T_ref_c>=T0):
		Cp_water = 4185.8518 * (0.9979 + 3.1e-6 * (T_ref_c - T0 - 35) ** 2 + 3.8e-9 * (T_ref_c - T0 - 35) ** 4)
	else:
		Cp_water = 4185.8518 * (1.000938 - 2.7052e-3 * (T_ref_c - T0) - 2.3235e-5 * (T_ref_c - T0) ** 2.	+ 4.3778e-6 * (T_ref_c - T0) ** 3 + 2.7136e-7 * (T_ref_c - T0) ** 4)
	e_sat_water_c=610.70 * m.exp(17.15 * (T_s_c - T0) / (T_s_c - 38.25))
	Le_c=4185.8518 * (597.3 - 0.561 * (T_s_c - T0))
	m_evap_per_c=0.622 * h_air[0, i_c]/ Cp_air * (e_sat_water_c/ (pressure[0, i_c] - e_sat_water_c) - e_inf / (p_0 - e_inf))
	m_evap_c=m_evap_per_c*area[0, i_c]
	m_imp_c = m_imp[n_c, i_c]
	m_in_c=0
	if(m_in[0, i_c-1]+m_imp_c<=m_evap_c):
		m_in_c=0
		m_evap_c=m_in_c+m_imp_c
		m_evap_per_c=m_evap_c/area[n_c, i_c]
	elif(T_wall_c + T0>=368.15):
		m_in_c=0
		m_evap_c=m_in[n_c, i_c-1]+m_imp_c
		m_evap_per_c=m_evap_c/area[n_c, i_c]
	else:
		m_in_c=(m_in[0, i_c-1]+m_imp_c-m_evap_c)
	m_out[n_c, i_c]=m_in_c
	Q_imp_per_c = m_imp_c * (u_inf * u_inf / 2 + Cp_water * T_inf)
	Q_imp_c = Q_imp_per_c * area[n_c, i_c]
	Q_evap_per_c = m_evap_per_c * (Le_c + Cp_water * T_s_c)
	Q_evap_c = Q_evap_per_c * area[n_c, i_c]
	Q_conv_per_c = h_air[n_c, i_c] * (T_s_c - 273.05)
	Q_conv_c=Q_conv_per_c*area[n_c, i_c]
	q_anti = 20150
	Q_anti_c= q_anti * area[n_c, i_c]

	Q_in_c=0

	T_ref_2 = (Q_in_c + Q_imp_c + Q_anti_c - Q_evap_c) / h_air[n_c, i_c] / area[n_c, i_c] + 247.8
	return T_ref_2





i=151

T_s_temp=303
T_ref_2_temp=289
T_wall_c=T_s_temp-T0
m_imp_per[0, i] = LWC * beta[0, i] * u_inf
m_imp[0, i] = m_imp_per[0, i] * area[0, i]
while (abs(T_s_temp-T_ref_2_temp)>0.00001):
	T_ref_2_temp=cell_iteration_c(151,0, T_wall_c, T_s_temp)
	T_s_temp=T_s_temp+0.01*(T_ref_2_temp-T_s_temp)

T_s=T_s_temp
print(T_s)









# def film_iteration( q_A,q_B, q_c, q_d,  NIND, c_start, c_end):
# 	n=c_end-c_start
# 	T_s = np.zeros((1, n))
# 	T_wall = np.zeros((1, n))
#
# 	T_s = np.zeros((1, n))
#
# 	Le = np.zeros((1, n))
# 	e_sat_water = np.zeros((1, n))
# 	m_evap_per = np.zeros((1, n))
# 	m_evap = np.zeros((1, n))
# 	m_in = np.zeros((1, n))
# 	m_out = np.zeros((1, n))
# 	m_stat_p = np.zeros((1, n))
# 	m_imp = np.zeros((1, n))
# 	m_imp_per = np.zeros((1, n))
# 	delt_w = np.zeros((1, n))
# 	Q_in = np.zeros((1, n))
# 	Q_out = np.zeros((1, n))
# 	Q_imp = np.zeros((1, n))
# 	Q_imp_per = np.zeros((1, n))
# 	Q_evap = np.zeros((1, n))
# 	Q_evap_per = np.zeros((1, n))
# 	Q_anti = np.zeros((1, n))
# 	Q_conv = np.zeros((1, n))
# 	Q_conv_per = np.zeros((1, n))
# 	Q_sta_p = np.zeros((1, n))
# 	q_anti = np.zeros((1, n))
#
#
# 	for i in range(c_start, c_end):
# 		if i ==c_start:
# 			T_s[0, i]=273.15