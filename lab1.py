import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.misc import derivative

n = 6
h = 1e-8
N = 10000

class Solution:
	def __init__(self, phis, vec):
		self.phis = phis
		self.vec = vec

	def __call__(self, x):
		res = 0
		for i in range(len(self.phis)):
			res += self.phis[i](x)*self.vec[i]
		return res

def integration(f1, f2, a, b, N=N):
	res = 0
	x = np.linspace(a, b, N)
	h = x[1] - x[0]
	for i in range(N):
		res += f1(x[i])*f2(x[i])*h
	return res

def df(x, func):
    return (func(x+h) - func(x-h))/(2*h)

def d2f(x, func):
	return (func(x+h) - 2*func(x) + func(x-h))/(h*h)

a,b,b1,b2,b3,c1,c2,c3,d1,d2,d3,k1,k2,p1,p2,q1,q2 \
 = 0,3,3,1,2,2,1,1,3,4,1,0,2,1,3,3,1
a1,a2,a3,a4,n1,n2,n3 = 1,3,2,5,7,5,2

def operatorA(f):
	k_df = lambda x: k(x)*df(x, f)
	# res = lambda x: -df(x, k)*df(x, f) - d2f(x, f)*k(x) + p(x)*df(x, f) + q(x)*f(x)
	res = lambda x: -k_df(x) + p(x)*df(x, f) + q(x)*f(x)

	return res 

def scalar_integral(f1, f2, a=a, b=b):
	func = lambda x: f1(x)*f2(x)
	return quad(func, a, b)[0]
	# return integration(f1, f2, a, b)

k = lambda x: b1*x**k1 + b2*x**k2 + b3
p = lambda x: c1*x**p1 + c2*x**p2 + c3
q = lambda x: d1*x**q1 + d2*x**q2 + d3

alpha = a1*a**n1 + a2*a**n2 + a3*a**n3 + a4
beta = -(a1*n1*a**(n1-1) + a2*n2*a**(n2-1) + \
	a3*n3*a**(n3-1))
gamma = a1*b**n1 + a2*b**n2 + a3*b**n3 + a4
delta = -(a1*n1*b**(n1-1) + a2*n2*b**(n2-1) + \
	a3*n3*b**(n3-1))

f = lambda x: -k(x)*(a1*n1*(n1-1)*x**(n1-2) + \
	a2*n2*(n2-1)*x**(n2-2) + a3*n3*(n3-1)*x**(n3-2)) - \
	(b1*k1*x**(k1-1)+b2*k2*x**(k2-1))*(a1*n1*x**(n1-1) +\
	a2*n2*x**(n3-1)+a3*n3*x**(n3-1)) + \
	p(x)*(a1*n1*x**(n1-1)+a2*n2*x**(n2-1) + \
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
var2 = alpha*(2*(x - b)*(x - B) + (x - b)**2) - beta*(phi2(x))
print(var1, var2)

phis = []
phis.append(phi1)
phis.append(phi2)
func = lambda i: lambda x: (x-b)*(x-b)*(x-a)**i

for i in range(2, n):
	phis.append(func(i))

Matr = np.zeros((n,n))
vec = np.array([.0 for i in range(n)])

for i in range(n):
	for j in range(n):
		Matr[i, j] = scalar_integral(operatorA(phis[j]), phis[i])
	vec[i] = scalar_integral(f, phis[i])


um_1 = lambda x, f: alpha*df(x, f) - beta*f(x)
um_2 = lambda x, f: gamma*df(x, f) + delta*f(x) 

sol = np.linalg.solve(Matr, vec)
print(Matr.dot(sol) - vec)
Sol = Solution(phis, sol)

print(um_1(a, Sol), um_2(b, Sol))

print(a, b, Sol(a), u_real(a), Sol(b), u_real(b))
x = np.linspace(a, b, N)
y = Sol(x)
y_real = u_real(x)
#print(y - y_real)
plt.plot(x, y, label='123')
plt.plot(x, y_real, label='real')
plt.legend()
plt.show()