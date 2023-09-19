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
abs_R=11
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


def R(abs_R,Bo):
	return 0.03895*abs_R*sqrt(Bo)

def E_abs(E,B):
	return 0.116*E*B/0.067



#def drexgap_asy(V, R, B, b):
#	sigma = 0.026185*(b**2)*(R**2)*V

#	return ( (0.116*B/b) + ( (8.9002*b*V/(sigma))*((1-(1/sqrt(sigma)))**2) ) )


B_range=np.arange(1.0,13.0,0.1)


E_ab_lis = map(lambda H: (E_abs(zero(0.067*V0/(0.116*H),R(abs_R,H),b,u,1.0),H)) - (E_abs(zero(0.067*V0/(0.116*H),R(abs_R,H),b,u,0.0),H)) , B_range)
#E_ab_lis_2 = map(lambda H: (E_abs(zero(0.067*V0/(0.116*H),R(abs_R,H),1.0,u,1.0),H)) - (E_abs(zero(0.067*V0/(0.116*H),R(abs_R,H),1.0,u,0.0),H)) , B_range)

# drexgap_asy_lis = map(lambda H: drexgap_asy(V0, abs_R, H, b) , B_range)


plt.plot(B_range, list(E_ab_lis), 'k-', label='Our model')
#plt.plot(B_range, list(E_ab_lis_2), 'k-', label=r'$\beta=1.0$')

# plt.plot(B_range, list(drexgap_asy_lis), 'k--', label=r'$(E_2-E_0)$, asymptotic')



x1=[float(j) for j in range(1,13)]

y1=[40.96774193548387,
	41.29032258064517,
	41.29032258064517,
	40.96774193548387,
	44.11290322580646,
	45.32258064516129,
	45.00000000000001,
	47.17741935483872,
	47.82258064516129,
	50.96774193548388,
	50.64516129032259,
	51.85483870967742]

plt.plot(x1,y1,'bo',label='Drexler et. al.')

plt.xlabel('Magnetic field B (T)')
plt.ylabel(r'$(E_2-E_0)$ (meV)')
#plt.title('(E_2 - E_0) vs B (B_i=B_o=B), R=8.85nm(best fit), beta=0.07, V=0.72eV')
plt.legend()
plt.show()
