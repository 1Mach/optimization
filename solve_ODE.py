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

def series(n):
	series_value=0
	for i in range(1, n):
		term_value=4*m.sin((2*n-1)*m.pi/2)  /  ((2*n-1)*m.pi)
		series_value+=term_value
	return series_value







def runge_kutta(y, x, dx, f):
	""" y is the initial value for y
        x is the initial value for x
        dx is the time step in x
        f is derivative of function y(t)
    """
	k1 = dx * f(y, x)
	k2 = dx * f(y + 0.5 * k1, x + 0.5 * dx)
	k3 = dx * f(y + 0.5 * k2, x + 0.5 * dx)
	k4 = dx * f(y + k3, x + dx)
	return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6.

if __name__=='__main__':
	t = 100
	x = 0.  #冰厚度5mm
	dt = 1
	xs, ts = [], []


def func(x, t):
	n=10
	series_value=0
	for i in range(1, n):
		term_value = 4 * m.sin((2 * n - 1) * m.pi / 2) / ((2 * n - 1) * m.pi)*m.exp((-alpha * (2*n-1)**2 * m.pi**2 * t)/(x**2))
		series_value += term_value
	return (1-series_value)*(h*(T_F-T_freestream)/k_l - P/(area*k_l))*(k_l/(rho_l*latent_heat_fusion))



while t <= 350 and t>=100:
	x = runge_kutta(x, t, dt, func)
	t += dt
	xs.append(x)
	ts.append(t)

exact = [(t ** 2 + 4) ** 2 / 16. for t in ts]
plt.plot(ts, xs, label='runge_kutta')
#plt.plot(ts, exact, label='exact')
plt.legend()
plt.show()