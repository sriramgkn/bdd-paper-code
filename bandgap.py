import scipy.special as sp
from math import*
import matplotlib.pyplot as plt
import numpy as np
import csv


V0=100.0 #100 meV
V=862.9*0.067 # 100 meV scaled wrt 0.116/0.067
u=1.0
dE=0.05
#abs_R=6
#b=1.0
#R=0.8
#l=0.0

def R(abs_R,B):
	return 0.03895*abs_R*sqrt(B)

def k1(E,b,u,l):
	return sqrt(2.0*b*E-u*(2.0*l+1.0))

def k2(E,V,R,b,u,l):
	return sqrt(-2*(E-V-l-0.5-(R**2.0*(u-1.0)/2.0)))

def f(E,V,R,b,u,l):
	ki=k1(E,b,u,l)

	ko=k2(E,V,R,b,u,l)

	x = ki*R*(sp.jvp(l,ki*R,1)/sp.jv(l,ki*R))

	return (b/2.0)+(b*ko*R)+(R**2*(b-u)/2.0)+x


def zero(V,R,b,u,l):
	E1=(l+0.5)*(u/b) + dE
	E=E1
	while f(E,V,R,b,u,l)*f(E+dE,V,R,b,u,l)>=0.0:
		E+=dE
	return (2.0*E+dE)/2.0

def asy_gap_01(R,V,b):
	rt_sigma = 1 + (b/2.0) + b*R*sqrt(0.067*0.026185*V) + 0.5*(b-1)*(R/25.7)**2
	return ( (38.19/(0.067*b*R**2))*(8.9)*(1.0-(1.0/rt_sigma))**2 ) - (0.116/(0.067*b))

def asy_gap_12(b):
	return (2*0.116/(0.067*b))

R_list=np.arange(0.35,0.8,0.005)

x=[r*25.67 for r in R_list]

#b_list = np.arange(0.1, 1.0, 0.01)

#print( 0.116*( zero(0.067*V0/(0.116),R(abs_R,1),0.7,1.0,-1.0) - zero(0.067*V0/(0.116),R(abs_R,1),0.7,1.0,0) )/0.067 )


# bandgap_data = list(map(lambda r: 0.116*(zero(V,r,1,u,-1.0) - zero(V,r,1,u,0.0)), R_list))
#
# with open('text.csv', 'w') as f:
# 	writer = csv.writer(f, delimiter='\t')
# 	writer.writerows(zip(x,bandgap_data))


plt.plot(x,list(map(lambda r: 0.116*(zero(V,r,0.7,u,-1.0) - zero(V,r,0.7,u,0.0))/0.067, R_list)),'b', label=r'$\beta = 0.7$, numerical')

plt.plot(x,list(map(lambda r: asy_gap_01(r,V0,0.7), x)),'b--', label=r'$\beta = 0.7 $, approx')

plt.plot(x,list(map(lambda r: 0.116*(zero(V,r,1.0,u,-1.0) - zero(V,r,1.0,u,0.0))/0.067, R_list)),'g', label=r'$\beta = 1.0$, numerical')

plt.plot(x,list(map(lambda r: asy_gap_01(r,V0,1.0), x)),'g--', label=r'$\beta = 1.0 $, approx')



plt.xlabel('Size R (nm)')
plt.ylabel(r'$(E_{1}-E_{0}) $ (meV)')
plt.legend(loc='best')

inset=plt.axes([0.62,0.38,0.25,0.25])


#for K in [8629.0, 43145.0]:
#	y_t = map(lambda b: 0.116*(zero(K,0.2,b,u,-1.0) - zero(K,0.2,b,u,0.0)), b_list)
#	y = list(y_t)
#	v_t = int(K/8629.0)
#	plt.plot(b_list, y, linewidth = 0.8, label=r'$V_{0} = %s $ eV' %v_t)

#plt.xlabel(r'$\beta$')
#plt.ylabel(r'$(E_{1}-E_{0}) $ (meV)')
#plt.legend(loc='best')

for b in [0.7,1.0]:
	y_t = map(lambda r: 0.116*(zero(V,r,b,u,1.0) - zero(V,r,b,u,-1.0))/0.067, R_list)
	y = list(y_t)
	if b==0.7:
		plt.plot(x,y,'b', label=r'$\beta = 0.7 $')
	else:
		plt.plot(x,y,'g', label=r'$\beta = 1.0 $')
plt.xlabel('Size R (nm)')
plt.ylabel(r'$(E_{2}-E_{1}) $ (meV)')
plt.legend(loc='center')

plt.show()
