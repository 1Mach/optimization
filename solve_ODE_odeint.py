import numpy as np
import sympy
import math as m
import matplotlib.pyplot as plt
from scipy import integrate
k_s=2.2
k_l=0.56
rho_s=917
rho_l=1000
c_ps=2010
alpha=k_s/c_ps  #thermal diffusivity
latent_heat_fusion=334000  #latent heat of fusion
h=34   #W/m^2-K
T_F=0 #fusion temperature
T_freestream=-5. #来流温度
P=1.64     #applied power(W) at -5 Celsius
area=0.3


def F(t, x):
	F=np.zeros(1)
	series_value=0
	n = 500
	for i in range(1, n):
		term_value = 4 * m.sin((2 * n - 1) * m.pi / 2) / ((2 * n - 1) * m.pi) * m.exp((-alpha * (2 * n - 1) ** 2 * m.pi ** 2 * t) / (x ** 2))  #发现除零了
		series_value += term_value
	F[0]=(1 - series_value) * (h * (T_F - T_freestream) / k_l - P / (area * k_l)) * (k_l / (rho_l * latent_heat_fusion))
	return F


def integrate(F, x, y, xStop, h):#四阶龙格库塔方法
	'''X,Y = integrate(F,x,y,xStop,h).
	4th-order Runge-Kutta method for solving the
	initial value problem {y}’ = {F(x,{y})}, where
	{y} = {y[0],y[1],...y[n-1]}.
	x,y = initial conditions
	xStop = terminal value of x
	h = increment of x used in integration
	F = user-supplied function that returns the
	array F(x,y) = {y’[0],y’[1],...,y’[n-1]}.
	'''
	def run_kut4(F, x, y, h):
		K0=h*F(x, y)
		K1=h*F(x+h/2.0, y+K0/2.0)
		K2=h*F(x+h/2.0, y+K1/2.0)
		K3=h*F(x+h, y+K2)
		return (K0+2.0*K1+2.0*K2+K3)/6.0
	X=[]
	Y=[]
	X.append(x)
	Y.append(y)
	while x<xStop:
		h=min(h, xStop-x)
		y=y+run_kut4(F, x, y, h)
		x=x+h
		X.append(x)
		Y.append(y)
	return np.array(X), np.array(Y)
t=100.0
tStop=350
x=np.array([0.])
err=0.1
freq=20
X,Y=integrate(F, t, x, tStop, err)

plt.plot(X, Y[:,0], 'o')
plt.show()