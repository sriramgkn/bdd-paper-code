import scipy.special as sp
from math import*
import matplotlib.pyplot as plt
import numpy as np
#x=np.arange(0.0,10.0,0.01)
#x=[j*(pi/100.0) for j in range(101)]
#x0=[sp.jv(0,j) for j in x]
#plt.plot(x,x0, label='initial state')
#plt.show()

V0=100.0
u=1.0
b=0.7
abs_R=11.0
#l=0.0

def k1(E,b,u,l):
	return sqrt(2.0*b*E-u*(2.0*l+1.0))

def k2(E,V,R,b,u,l):
	return sqrt(-2*(E-V-l-0.5-(R**2.0*(u-1.0)/2.0)))

def f(E,V,R,b,u,l):
	ki=k1(E,b,u,l)

	ko=k2(E,V,R,b,u,l)

	x = ki*R*(sp.jvp(l,ki*R,1)/sp.jv(l,ki*R))

	return (b/2.0)+(b*ko*R)+(R**2*(b-u)/2.0)+x

dE=0.01

def zero(V,R,b,u,l):
	E1=(l+0.5)*(u/b)+dE
	E=E1
	while f(E,V,R,b,u,l)*f(E+dE,V,R,b,u,l)>=0.0:
		E+=dE
	return (2.0*E+dE)/2.0


def R(abs_R,B):
	return 0.03895*abs_R*sqrt(B)

def E_abs(E,B):
	return 0.116*E*B/0.067

#Asymptotic functions

def asy_E_0(R,V,B,b):
	rt_sigma = 1 + (b/2.0) + b*R*sqrt(0.067*0.026185*V) + 0.5*B*(b-1)*(R/25.7)**2
	return ( (38.19/(0.067*b*R**2))*(2.405**2)*(1.0-(1.0/rt_sigma))**2 ) + (0.5*0.116*B/(0.067*b))

def asy_E_1(R,V,B,b):
	rt_sigma = 1 + (b/2.0) + b*R*sqrt(0.067*0.026185*V) + 0.5*B*(b-1)*(R/25.7)**2
	return ( (38.19/(0.067*b*R**2))*(3.832**2)*(1.0-(1.0/rt_sigma))**2 ) + (1.5*0.116*B/(0.067*b))

def asy_E_11(R,V,B,b):
	rt_sigma = 1 + (b/2.0) + b*R*sqrt(0.067*0.026185*V) + 0.5*B*(b-1)*(R/25.7)**2
	return ( (38.19/(0.067*b*R**2))*(3.832**2)*(1.0-(1.0/rt_sigma))**2 ) - (0.5*0.116*B/(0.067*b))

def asy_E_2(R,V,B,b):
	rt_sigma = 1 + (b/2.0) + b*R*sqrt(0.067*0.026185*V) + 0.5*B*(b-1)*(R/25.7)**2
	return ( (38.19/(0.067*b*R**2))*(5.136**2)*(1.0-(1.0/rt_sigma))**2 ) + (2.5*0.116*B/(0.067*b))

def asy_E_22(R,V,B,b):
	rt_sigma = 1 + (b/2.0) + b*R*sqrt(0.067*0.026185*V) + 0.5*B*(b-1)*(R/25.7)**2
	return ( (38.19/(0.067*b*R**2))*(5.136**2)*(1.0-(1.0/rt_sigma))**2 ) - (1.5*0.116*B/(0.067*b))


B_range = np.arange(0.1,10,0.2)

#fig, ax1 = plt.subplots()
#x_cor, y_cor, width, height = [0.17, 0.67, 0.2, 0.2]
#ax2 = fig.add_axes([x_cor, y_cor, width, height])

for l in [0, 1, -1, 2, -2]:
	E_ab_list_t = map(lambda H: E_abs(zero(0.067*V0/(0.116*H),R(abs_R,H),b,u,l),H) , B_range)
	E_ab_list = list(E_ab_list_t)
	plt.plot(B_range, E_ab_list, label='(1,%s), numerical' % l)


E_as_0 = map(lambda H: asy_E_0(abs_R,V0,H,b) , B_range)
E_as_1 = map(lambda H: asy_E_1(abs_R,V0,H,b) , B_range)
E_as_11 = map(lambda H: asy_E_11(abs_R,V0,H,b) , B_range)
E_as_2 = map(lambda H: asy_E_2(abs_R,V0,H,b) , B_range)
E_as_22 = map(lambda H: asy_E_22(abs_R,V0,H,b) , B_range)

plt.gca().set_prop_cycle(None)

plt.plot(B_range, list(E_as_0),'--', label='(1,0), approx') # #1f77b4
plt.plot(B_range, list(E_as_1),'--', label='(1,1), approx') # #ff7f0e
plt.plot(B_range, list(E_as_11),'--', label='(1,-1), approx') # #2ca02c
plt.plot(B_range, list(E_as_2),'--', label='(1,2), approx') # #d62728
plt.plot(B_range, list(E_as_22),'--', label='(1,-2), approx') # #9467bd

plt.xlabel('Magnetic field B (T)')
plt.ylabel('Energy (meV)')
plt.legend(loc=1)

inset=plt.axes([0.17, 0.67, 0.2, 0.2])

for l in [0, 1, -1, 2, -2]:
	E_ab_list_t = map(lambda H: E_abs(zero(V0/(0.116*H),R(abs_R,H),1.0,u,l),H) , B_range)
	E_ab_list = list(E_ab_list_t)
	plt.plot(B_range, E_ab_list)

plt.text(5,160,r'$\beta=0.7$')
plt.show()
