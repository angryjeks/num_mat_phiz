import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

n = 5

a,b,b1,b2,b3,c1,c2,c3,d1,d2,d3,k1,k2,p1,p2,q1,q2 \
 = 0,3,3,1,2,2,1,1,3,4,1,0,2,1,3,3,1
a1,a2,a3,a4,n1,n2,n3 = 1,3,2,5,7,5,2

k = lambda x: b1*x**k1 + b2*x**k2 + b3
p = lambda x: c1*x**p1 + c2*x**p2 + c3
q = lambda x: d1*x**q1 + d2*x**q2 + d3

alpha = a1*a**n1 + a2*a**n2 + a3*a**n3 + a4
beta = -(a1*n1*a**(n1-1) + a2*n2*a**(n2-1) + 
	a3*n3*a**(n3-1))
gamma = a1*b**n1 + a2*b**n2 + a3*b**n3 + a4
delta = -(a1*n1*b**(n1-1) + a2*n2*b**(n2-1) + 
	a3*n3*b**(n3-1))

f = lambda x: -k(x)*(a1*n1*(n1-1)*x**(n1-2) + \
	a2*n2*(n2-1)*x**(n2-2) + a3*n3*(n3-1)*x**(n3-2)) - \
	(b1*k1*x**(k1-1)+b2*k2*x**(k2-1))*(a1*n1*x**(n1-1) +\
	a2*n2*x**(n3-1)+a3*n3*x**(n3-1)) + \
	p(x)*(a1*n1*x**(n1-1)+a2*n2*x**(n2-1)+ \
	a3*n3*x**(n3-1)) + q(x)*(a1*x**n1+a2*x**n2 + \
	a3*x**n3+a4)

u_real = lambda x: a1*x**n1 + a2*x**n2 + a3*x**n3 + a4

A = gamma*(b-a)/(2*gamma + delta*(b-a)) + b
B = a - alpha*(a-b)/(beta*(a-b)-2*alpha)

phi1 = lambda x: (x-a)*(x-a)*(x-A)
phi2 = lambda x: (x-b)*(x-b)*(x-B)
x = b
var1 = gamma*(2*(x - a)*(x - A) + (x - a)**2) + delta*(phi1(x))
x = a
var2 = alpha*(2*(x - b)*(x - B) + (x - b)**2) + beta*(phi2(x))
print(var1, var2)

phis = []
phis.append(phi1)
phis.append(phi2)
func = lambda i: lambda x: (x-b)*(x-b)*(x-a)**i

for i in range(2, n):
	phis.append(func(i))

#print(phis[2](2), phis[3](2))









#x = np.linspace(0, 2, 1000)
#y = u_real(x)
#plt.plot(x,y)
#plt.show()