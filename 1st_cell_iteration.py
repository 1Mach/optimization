import numpy as np
import math as m



n=6 #网格个数


h_s=np.zeros((n, 1))
h_s[0,0]=296.4705882
h_s[1, 0]=248.0672269
h_s[2,0]=235.9663866
h_s[3, 0]=222.3529412
h_s[4,0]=211.7647059
h_s[5, 0]=193.6134454

p=np.zeros((n, 1))
p[0, 0]=74959.70313
p[1, 0]=74949.42969
p[2, 0]=74923.97656
p[3, 0]=74883.60156
p[4, 0]=74826.33594
p[5, 0]=74753.32031




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


#-------------------液膜表面饱和水蒸气压---------------------
e_inf=610.70*m.exp(17.15*(T_inf-T0)/(T_inf-38.25))






T_s=np.zeros((n, 1))  #
T_wall=np.zeros((n, 1))
T_rec=np.zeros((n, 1))
e_sat_water=np.zeros((n, 1))
Le=np.zeros((n, 1))
m_evap_per=np.zeros((n, 1))
m_evap=np.zeros((n, 1))
m_in=np.zeros((n, 1))
m_out=np.zeros((n, 1))
m_stat_p=np.zeros((n, 1))




i=0
T_s[i, 0]=273.15
T_wall[i, 0]=T_s[i, 0]-T0
#T_rec[i, 0]=T_inf+(u_inf**2)/2/Cp_air*(1-((ue(i)/u_inf)^2)*(1-0.85))

for j in range(0, 1000):
	T_ref=0.5*(T_s[i, 0]+T_wall[i, 0]+T0)
	if(T_ref>=T0):
		row_water=(999.8396+18.224944*(T_ref-T0)-7.922210e-3*(T_ref-T0)**2.-55.44846e-6*(T_ref-T0)^3+149.7562e-9*(T_ref-T0)^4-393.2952e-12*(T_ref-T0)^5)/(1+(18.159725e-3)*(T_ref-T0))
		Cp_water=4185.8518*(0.9979+3.1e-6*(T_ref-T0-35)**2+3.8e-9*(T_ref-T0-35)**4)
		miu_water=1e-3*(1.76*m.exp(-3.5254e-2*(T_ref-T0)+4.7163e-4*(T_ref-T0)**2.-6.0667e-6*(T_ref-T0)**3))
	else:
		rou_water = 1000*(0.99986+6.69e-5*(T_ref-T0)-8.486e-6*(T_ref - T0) ** 2.	+ 1.518e-7 * (T_ref - T0) ** 3 - 6.9484e-9 * (T_ref - T0) ** 4 - 3.6449e-10 ** (T_ref - T0) ** 5.- 7.497e-12 * (T_ref -T0)**6 T0) ** 6)
		Cp_water = 4185.8518 * (1.000938 - 2.7052e-3 * (T_ref - T0) - 2.3235e-5 * (T_ref - T0) ** 2.	+ 4.3778e-6 * (T_ref - T0) ** 3 + 2.7136e-7 * (T_ref - T0) ** 4)
		miu_water = 1e-3 * (1.76 * m.exp(-5.5721e-2 * (T_ref - T0) - 1.3943e-3 * (T_ref - T0) ** 2. - 4.3015e-5 * (T_ref - T0) ** 3))
	#饱和蒸汽压
	e_sat_water[i, 0] = 610.70 * m.exp(17.15 * (T_s[i, 0] - T0) / (T_s[i, 0] - 38.25))

	#蒸发热
	Le[i, 0]= 4185.8518 * (597.3 - 0.561 * (T_s[i, 0] - T0))
	m_evap_per[i, 0] = 0.622 * h(i) / Cp_air * (e_sat_water[i, 0] / (p[i, 0] - e_sat_water[i, 0]) - e_inf / (p_0 - e_inf))
	m_evap[i, 0]=m_evap_per[i, 0]*delt_s  #这个需要替换成我自己的结果
	m_in[i, 0]=0
	#下面需要继续编写





