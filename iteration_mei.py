import pandas as pd
import numpy as np
from scipy.interpolate import pchip_interpolate
import math as m
import matplotlib.pyplot as plt
#把所有量都换成行的形式

#read the original results as a two-dimentional array like
# [[a]
#  [b]
#  [c]]
s=np.transpose(np.asarray(pd.read_csv("./original_data/fixed Twall-SA.csv",skiprows=[0], usecols=[2] )))
#print(s)

p_s=np.transpose(np.asarray(pd.read_csv("./original_data/fixed Twall-SA.csv",skiprows=[0], usecols=[3] )))
p_coef_s=np.transpose(np.asarray(pd.read_csv("./original_data/fixed Twall-SA.csv",skiprows=[0], usecols=[4] )))
beta_s=np.transpose(np.asarray(pd.read_csv("./original_data/fixed Twall-SA.csv",skiprows=[0], usecols=[5] )))
F_s=np.transpose(np.asarray(pd.read_csv("./original_data/fixed Twall-SA.csv",skiprows=[0], usecols=[6] )))
ue_s=np.transpose(np.asarray(pd.read_csv("./original_data/fixed Twall-SA.csv",skiprows=[0], usecols=[8] )))
s_c_s=np.transpose(np.asarray(pd.read_csv("./original_data/h-silva.csv",skiprows=[0], usecols=[0] ).dropna()))
h_s=np.transpose(np.asarray(pd.read_csv("./original_data/h-silva.csv",skiprows=[0], usecols=[1], skip_blank_lines=True).dropna()))  #drop NaN Value


s_exp=np.transpose(np.asarray(pd.read_csv("./original_data/Twall-exp.csv",skiprows=[0], usecols=[0], skip_blank_lines=True).dropna()))
t_exp=np.transpose(np.asarray(pd.read_csv("./original_data/Twall-exp.csv",skiprows=[0], usecols=[1], skip_blank_lines=True).dropna()))

s_antice=np.transpose(np.asarray(pd.read_csv("./original_data/Twall-antice.csv",skiprows=[0], usecols=[0], skip_blank_lines=True).dropna()))
t_antice=np.transpose(np.asarray(pd.read_csv("./original_data/Twall-antice.csv",skiprows=[0], usecols=[1], skip_blank_lines=True).dropna()))
#------------------------------------------------------------------------------------------------------------------


#----------------------fitting the data to s calculated by CFD---------------------------------------------------
x_s=np.arange(-0.18, 0.18, 1e-4)

xx_s=x_s.reshape((1, -1))#-1表示任意行数，1表示1列，相当于矩阵转置
#print(xx_s)
yy_p_coef_s=(pchip_interpolate(s[0,:], p_coef_s[0,:], xx_s[0,:] )).reshape((1, -1))
#plt.plot(xx_s[0, :], yy_p_coef_s[0, :])
yy_p_s=(pchip_interpolate(s[0,:], p_s[0,:], xx_s[0,:] )).reshape((1, -1))
yy_beta_s=(pchip_interpolate(s[0,:], beta_s[0,:], xx_s[0,:] )).reshape((1, -1))
yy_F_s=(pchip_interpolate(s[0,:], F_s[0,:], xx_s[0,:] )).reshape((1, -1))
yy_ue_s=(pchip_interpolate(s[0,:], ue_s[0,:], xx_s[0,:] )).reshape((1, -1))
n_s=np.size(xx_s, 1)#计算xx_s第1轴（多少列）的元素个数
#print(n_s)
A_s=np.vstack((xx_s, yy_p_s, yy_p_coef_s, yy_beta_s, yy_F_s, yy_ue_s))  #每一列代表一个参数
#print(A_s)
num=(np.argmin(yy_F_s, axis=1))[0]
#print(num)
a_s=(A_s[ :,0:num])
b_s=a_s[:, ::-1]  #将列倒过来排列；下表面
#print(b_s)
c_s=A_s[:,num:n_s]
#print(c_s)
s_up=c_s[[0], :]  #矩阵
#print(s_up)
p_up=c_s[[1], :]
p_coef_up=c_s[[2], :]
beta_up=c_s[[3], :]
F_up=c_s[[4], :]
ue_up=c_s[[5], :]

s_down=b_s[[0],: ]
p_down=b_s[[1],: ]
p_coef_down=b_s[[2],: ]
beta_down=b_s[[3],: ]
F_down=b_s[[4],: ]
ue_down=b_s[[5],: ]

x_s0=np.arange(-0.18, 0.18, 1e-4)
xx_s0=x_s.reshape((1, -1))
s_0=s_c_s * 0.914
yy_h_s=pchip_interpolate(s_0[0, :], h_s[0, :], xx_s0[0, :])
n_s0=np.size(xx_s0, 1)
A_s0=np.vstack((xx_s0, yy_h_s))
s_cri=s_up[0, 0]  #中间点

for i_0 in range(0, n_s0):
	if (xx_s0[0, i_0]>=s_cri):
		break
	
a_s0=A_s0[:,0:i_0]  #下表面
b_s0=a_s0[:, ::-1 ]  #上表面
c_s0=A_s0[:,i_0:n_s0]

s0_up=c_s0[[0],:]
h_up=c_s0[[1],:]
s0_down=b_s0[[0],:]
h_down=b_s0[[1],:]

q_anti_s=np.zeros((1, n_s ))
S_C_s=xx_s/0.914

x_s_antice=np.arange(-0.137, 0.137, 1e-3)
xx_s_antice=x_s_antice.reshape((1, -1))
yy_s_antice=pchip_interpolate(s_antice[0, :], t_antice[0, :], xx_s_antice[0, :])
n_s_antice=np.size(xx_s_antice, 1)
A_s_antice=np.vstack((xx_s_antice, yy_s_antice))



for i_antice in range(0, n_s_antice):
	if (xx_s_antice[0,i_antice]>=s_cri):
		break
#print(i_antice)

a_antice=A_s_antice[ :,0:i_antice]
b_antice=a_antice[:,::-1]
c_antice=A_s_antice[:,i_antice: n_s_antice]

s_antice_up=c_antice[[0],:]
t_antice_up=c_antice[[1],:]
s_antice_down=b_antice[[0],:]
t_antice_down=b_antice[[1],:]



print(n_s)
for i in range(0, n_s):
	if ((xx_s[0, i]>=-0.093599) and (xx_s[0, i]<-0.055499)):
		q_anti_s[0, i]=20.15*1e3
	elif ((xx_s[0, i]>=-0.055499) and (xx_s[0, i]<-0.030099)):
		q_anti_s[0, i]=21.70*1e+3
	elif ((xx_s[0, i]>=-0.030099) and (xx_s[0, i]<-0.004699)):
		q_anti_s[0, i]=32.55*1e+3
	elif ((xx_s[0, i]>=-0.004699) and (xx_s[0, i]<0.014351)):
		q_anti_s[0, i]=43.4*1e+3
	elif ((xx_s[0, i]>=0.01435) and (xx_s[0, i]<0.039751)):
		q_anti_s[0, i]=26.35*1e+3
	elif ((xx_s[0, i]>=0.039751) and (xx_s[0, i]<0.065151)):
		q_anti_s[0, i]=18.60*1e+3
	elif ((xx_s[0, i]>=0.065151) and (xx_s[0, i]<0.103251)):
		q_anti_s[0, i]=18.60*1e+3
	else:
		q_anti_s[0, i]=0
#print(q_anti_s)
A_s1=np.vstack((xx_s, q_anti_s))

a_s1=A_s1[:, 0:num]
b_s1=a_s1[:, ::-1]
c_s1=A_s1[:, num:n_s]

q_anti_down=b_s1[[1],:]
q_anti_up=c_s1[[1],:]

f_s=S_C_s[:, 0:num]
S_C_down=f_s[:, ::-1]
S_C_up=S_C_s[:, num:n_s]

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

#-----------------液膜表面饱和水蒸气压------------
e_inf=610.70*m.exp(17.15*(T_inf-T0)/(T_inf-38.25)) #本来是(T_inf-T0(T0是冻结温度))

h=h_up
p=p_up
beta=beta_up
F_shear=F_up
q_anti=q_anti_up
ue=ue_up

n=np.size(h, 1)  #h列的个数，也就是网格的个数

#print(n)


e_sat_water=np.zeros((1, n))

Le=np.zeros((1, n))

m_evap_per=np.zeros((1, n))
m_evap=np.zeros((1, n))
m_in=np.zeros((1, n))
m_out=np.zeros((1, n))

delt_p=np.zeros((1, n))
delt_w=np.zeros((1, n))
us=np.zeros((1, n))
err=1e-20

Q_imp_per=np.zeros((1, n))
Q_imp=np.zeros((1, n))
Q_evap_per=np.zeros((1, n))
Q_evap=np.zeros((1, n))
Q_conv_per=np.zeros((1, n))
Q_conv=np.zeros((1, n))
Q_in=np.zeros((1, n))
Q_out=np.zeros((1, n))
Q_anti=np.zeros((1, n))
Q_wall=np.zeros((1, n))
q_wall=np.zeros((1, n))
T_s=np.zeros((1, n))
T_wall=np.zeros((1, n))
T_rec=np.zeros((1, n))
delt_T=np.zeros((1, n))
q_anti2=np.ones((1, n))
m_sta_p=np.ones((1, n))
Q_sta_p=np.ones((1, n))
#print(q_anti2)

m_imp_per=LWC * beta *u_inf
m_imp=m_imp_per * delt_s

num_dry=-1

for i in range(0, n):
	if i==0:
		T_s[0, i]=273.15
		T_wall[0, i]=T_s[0, i]- T0
		us[0, i]=0
		T_rec[0, i]=T_inf + (u_inf**2)/2/Cp_air*(1-((ue[0, i]/u_inf)**2)*(1-0.85))
		print(T_rec[0, i])
		for j in range(0, 1000):
			T_ref=0.5 * (T_s[0, i]+T_wall[0, i]+T0)
			if (T_ref>=T0):
				rou_water=(999.8396+18.224944*(T_ref-T0)-7.922210e-3*(T_ref-T0)**2.-55.44846e-6*(T_ref-T0)**3+149.7562e-9*(T_ref-T0)**4-393.2952e-12*(T_ref-T0)**5)/(1+(18.159725e-3)*(T_ref-T0))
				Cp_water=4185.8518*(0.9979+3.1e-6*(T_ref-T0-35)**2+3.8e-9*(T_ref-T0-35)**4)
				miu_water=1e-3*(1.76*m.exp(-3.5254e-2*(T_ref-T0)+4.7163e-4*(T_ref-T0)**2.-6.0667e-6*(T_ref-T0)**3))
			else:
				rou_water = 1000 * (0.99986 + 6.69e-5 * (T_ref - T0) - 8.486e-6 * (T_ref - T0) ** 2.	+ 1.518e-7 * (T_ref - T0) ** 3 - 6.9484e-9 * (T_ref - T0) ** 4 - 3.6449e-10 ** (T_ref - T0) ** 5.- 7.497e-12 * (T_ref - T0) ** 6)
				Cp_water = 4185.8518 * (1.000938 - 2.7052e-3 * (T_ref - T0) - 2.3235e-5 * (T_ref - T0) ** 2.	+ 4.3778e-6 * (T_ref - T0) ** 3 + 2.7136e-7 * (T_ref - T0) ** 4)
				miu_water = 1e-3 * (1.76 * m.exp(-5.5721e-2 * (T_ref - T0) - 1.3943e-3 * (T_ref - T0) ** 2. - 4.3015e-5 * (T_ref - T0) ** 3))
			D_v=0.211*((0.5*(T_ref+T_inf)/273.15)**1.94)*(101325/p_0)*1e-4
			Sc = miu_air/rou_air/D_v
			
			e_sat_water[0, i ] = 610.70 * m.exp(17.15 * (T_s[0, i] - T0) / (T_s[0, i] - 38.25))
			Le[0, i]= 4185.8518 * (597.3 - 0.561 * (T_s[0, i] - T0))
			m_evap_per[0, i]=0.622 * h[0, i]/ Cp_air * (e_sat_water[0, i] / (p[0, i] - e_sat_water[0, i]) - e_inf / (p_0 - e_inf))
			m_evap[0, i]=m_evap_per[0, i]*delt_s
			m_in[0, i]=0
			if (m_in[0, i]+m_imp[0, i]<=m_evap[0, i]):
				m_sta_p[0, i]=0
				m_evap[0, i]=m_in[0, i]+m_imp[0, i]
				m_evap_per[0, i]=m_evap[0, i]/delt_s
			elif(T_wall[0, i]+T0>=368.15):
				m_sta_p[0, i]=0
				m_evap[0, i]=m_in[0, i]+m_imp[0, i]
				m_evap_per[0, i] = m_evap[0, i] / delt_s
			else:
				m_sta_p[0, i]=m_imp[0, i]-m_evap[0, i]
			
			delt_p[0, i]=(p[0, i+1]-p[0, i])/delt_s
			if (m_sta_p[0, i]==0):
				delt_w[0, i]=0
			else:
				m_sta_p[0 ,i]=m_imp[0 ,i]-m_evap[0 ,i]
			if (m_sta_p[0 ,i]==0):
				delt_w[0, i]=0
			else:
				x0=1e-6
				for l in range(0, 1000):
					f=-2/3*delt_p[0, i]*x0**3-1/6*F_shear[0, i]/miu_water*m_imp_per[0, i]*x0**3+F_shear[0, i]*x0**2 +2/3*m_imp_per[0, i]*u_inf*x0**2-m_imp_per[0, i]*m_sta_p[0, i]/rou_water*x0-2*m_sta_p[0, i]*miu_water/rou_water
					df=-2*delt_p[0, i]*x0**2-1/2*F_shear[0, i]/miu_water*m_imp_per[0, i]*x0**2+2*F_shear[0, i]*x0 +4/3*m_imp_per[0, i]*u_inf*x0-m_imp_per[0, i]*m_sta_p[0, i]/rou_water
					x1=x0-f/df
					if (abs(x1-x0)<err):
						break
					x0=x1
				delt_w[0, i]=x0
			Q_in[0, i]=0
			Q_imp_per[0, i]=m_imp_per[0, i]*(u_inf*u_inf/2 + Cp_water*T_inf)
			Q_imp[0, i]=Q_imp_per[0, i]*delt_s
			Q_evap_per[0, i]=m_evap_per[0, i]*(Le[0, i] +Cp_water*T_s[0, i])
			Q_evap[0, i]=Q_evap_per[0, i]*delt_s
			Q_anti[0, i]=q_anti[0, i]*delt_s
			
			Q_conv_per[0, i]=h[0, i]*(T_s[0, i]-T_rec[0, i])
			Q_conv[0, i]=Q_conv_per[0, i]*delt_s
			Q_sta_p[0, i]=(Q_imp[0, i]+Q_anti[0, i])-(Q_evap[0, i]+Q_conv[0, i])
			
			T_ref_1=Q_sta_p[0, i]/m_sta_p[0, i]/Cp_water
			
			if ((abs(T_s[0, i]-T_ref_1))<=0.00001):
				break
			else:
				T_s[0, i]=T_s[0, i]+0.01*(T_ref_1-T_s[0, i])
			T_wall[0, i]=T_s[0, i]+q_anti[0, i]*delt_w[0, i]/lamda_water-T0
	else:
		if (num_dry==-1):
			T_s[0, i]=T_s[0, i-1]
			T_wall[0, i]=T_wall[0, i-1]
			T_rec[0, i]=T_inf + (u_inf**2)/2/Cp_air*(1-((ue[0, i]/u_inf)**2)*(1-0.89))
			for j_1 in range(0, 100000):
				T_ref=0.5 * (T_s[0, i]+T_wall[0, i]+T0)
				if (T_ref>=T0):
					rou_water=(999.8396+18.224944*(T_ref-T0)-7.922210e-3*(T_ref-T0)**2.-55.44846e-6*(T_ref-T0)**3+149.7562e-9*(T_ref-T0)**4-393.2952e-12*(T_ref-T0)**5)/(1+(18.159725e-3)*(T_ref-T0))
					Cp_water=4185.8518*(0.9979+3.1e-6*(T_ref-T0-35)**2+3.8e-9*(T_ref-T0-35)**4)
					miu_water=1e-3*(1.76*m.exp(-3.5254e-2*(T_ref-T0)+4.7163e-4*(T_ref-T0)**2.-6.0667e-6*(T_ref-T0)**3))
				else:
					rou_water = 1000 * (0.99986 + 6.69e-5 * (T_ref - T0) - 8.486e-6 * (T_ref - T0) ** 2.	+ 1.518e-7 * (T_ref - T0) ** 3 - 6.9484e-9 * (T_ref - T0) ** 4 - 3.6449e-10 * (T_ref - T0) ** 5.- 7.497e-12 * (T_ref - T0) ** 6)
					Cp_water = 4185.8518 * (1.000938 - 2.7052e-3 * (T_ref - T0) - 2.3235e-5 * (T_ref - T0) ** 2. + 4.3778e-6 * (T_ref - T0) ** 3 + 2.7136e-7 * (T_ref - T0) ** 4)
					miu_water = 1e-3 * (1.76 * m.exp(-5.5721e-2 * (T_ref - T0) - 1.3943e-3 * (T_ref - T0) ** 2. - 4.3015e-5 * (T_ref - T0) ** 3))
				D_v=0.211*((0.5*(T_ref+T_inf)/273.15)**1.94)*(101325/p_0)*1e-4
				Sc = miu_air/rou_air/D_v
				e_sat_water[0, i ] = 610.70 * m.exp(17.15 * (T_s[0, i] - T0) / (T_s[0, i] - 38.25))
				Le[0, i]= 4185.8518 * (597.3 - 0.561 * (T_s[0, i] - T0))
				m_evap_per[0, i]=0.622 * h[0, i]/ Cp_air * (e_sat_water[0, i] / (p[0, i] - e_sat_water[0, i]) - e_inf / (p_0 - e_inf))
				m_evap[0, i]=m_evap_per[0, i]*delt_s
				if (m_in[0, i-1]+m_imp[0, i]<=m_evap[0, i]):
					m_in[0, i]=0
					m_evap[0, i]=m_in[0, i-1]+m_imp[0, i]
					m_evap_per[0, i]=m_evap[0, i]/delt_s
				elif(T_wall[0, i]+T0>=368.15):
					m_in[0, i]=0
					m_evap[0, i]=m_in[0, i-1]+m_imp[0, i]
					m_evap_per[0, i] = m_evap[0, i] / delt_s
				else:
					m_in[0, i]=m_in[0, i-1]+m_imp[0, i]-m_evap[0, i]
				m_out[0, i-1]=m_in[0, i]
				
				delt_p[0, i]=(p[0, i+1]-p[0, i-1])/delt_s/2
				if(m_in[0, i]==0):
					delt_w[0, i]=0
					us[0, i]=0
					num_dry_1=i
					if (m_imp[0, i+1]>0):
						num_dry=-1  #撞击区内干防冰
					else:
						num_dry=0  #撞击区外
					Q_imp_per[0, i]=m_imp_per[0, i]*(u_inf*u_inf/2 + Cp_water*T_inf)
					Q_imp[0, i]=Q_imp_per[0, i]*delt_s
					Q_anti[0, i]=q_anti[0, i] * delt_s
					for num_j in range(1, 2):
						Q_evap_per[0, i]=m_evap_per[0, i]*(Le[0, i] +Cp_water*T_s[0, i])
						Q_evap[0, i]=Q_evap_per[0, i]*delt_s
						Q_conv_per[0, i]=h[0, i]*(T_s[0, i]-T_rec[0, i])
						Q_conv[0, i]=Q_conv_per[0, i]*delt_s
						T_ref_3=(Q_in[0, i-1]+Q_imp[0, i]+Q_anti[0, i] -Q_evap[0, i])/h[0, i]/delt_s + T_rec[0, i]
					T_wall[0, i]=T_ref_3-T0
					break
				else:
					x0=1e-6
					
					for l in range(0, 1000):
						f=-2/3*delt_p[0, i]*x0**3-1/6*F_shear[0, i]/miu_water*m_imp_per[0, i]*x0**3+F_shear[0, i]*x0**2.+2/3*m_imp_per[0, i]*u_inf*x0**2-m_imp_per[0, i]*m_in[0, i]/rou_water*x0-2*m_in[0, i]*miu_water/rou_water
						df=-2*delt_p[0, i]*x0**2-1/2*F_shear[0, i]/miu_water*m_imp_per[0, i]*x0**2+2*F_shear[0, i]*x0+4/3*m_imp_per[0, i]*u_inf*x0-m_imp_per[0, i]*m_in[0, i]/rou_water
						x1=x0-f/df
						if (abs(x1-x0)<err):
							break
						x0=x1
					delt_w[0, i]=x0
					us[0, i]=(-delt_w[0, i]*delt_p[0, i]+2*F_shear[0, i]+m_imp_per[0, i]*u_inf)/(2*miu_water/delt_w[0, i]+m_imp_per[0, i])
					
					Q_imp_per[0, i]=m_imp_per[0, i]*(u_inf*u_inf/2 + Cp_water*T_inf)
					Q_imp[0, i]=Q_imp_per[0, i]*delt_s
					Q_evap_per[0, i]=m_evap_per[0, i]*(Le[0, i] +Cp_water*T_s[0, i])
					Q_evap[0, i]=Q_evap_per[0, i]*delt_s
					Q_conv_per[0, i]=h[0, i]*(T_s[0, i]-T_rec[0, i])
					Q_conv[0, i]=Q_conv_per[0, i]*delt_s
					Q_anti[0, i]=q_anti[0, i]*delt_s
					
					Q_in[0, i]=Q_in[0, i-1]+(Q_imp[0, i]+Q_anti[0, i])-(Q_evap[0, i]+Q_conv[0, i])
					Q_out[0, i-1]=Q_in[0, i]
					T_ref_2=(Q_in[0, i]/rou_water/Cp_water-(1/4*us[0, i]/delt_w[0, i]-1/12*F_shear[0, i]/miu_water)*q_anti[0, i]/lamda_water*(delt_w[0, i]**3))/((2/3*us[0, i]/delt_w[0, i]-1/6*F_shear[0, i]/miu_water)*(delt_w[0, i]**2))
					if (abs(T_s[0, i]-T_ref_2)<=0.00001):
						break
					else:
						T_s[0, i]=T_s[0, i]+0.01 * (T_ref_2-T_s[0, i])
					T_wall[0, i]=T_s[0, i]+q_anti[0, i]*delt_w[0, i]/lamda_water-T0
		else:
			break
	delt_T[0, i]=T_wall[0, i]+T0-T_s[0, i]

T_wall_2=T_wall[0]
				
					
				
				
				
				
				
				
				
				
				
				
				
				
				
				
