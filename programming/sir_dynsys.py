import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

a = 2
#b = 1

def plotter(b):
	func = lambda x : a - b*x - np.exp(-x)
	guess = a/b
	max_x = fsolve(func, guess)
	x = np.arange(0.0, max_x*1.05, 0.01)
	y1 = a - b*x
	y2 = np.exp(-x)
	y3 = y1 - y2
	null = 0*x
	plt.figure()
	plt.fill_between(x, y1, y2)
	plt.plot(x, y3, color="yellow", label="difference")
	plt.plot(x, null, "--b")
	plt.title('b = ' + str(b))
	plt.legend(loc="best")
	plt.xlabel('time')
	plt.ylabel('P(t)')
	plt.savefig("fig/SIR_b_" + str(b) + ".pdf")

def stability_plotter(b):
	func = lambda x : a - b*x - np.exp(-x)
	guess = a/b
	max_x = fsolve(func, guess)
	x = np.arange(0.0, max_x*1.2, 0.01)
	y = a - b*x - np.exp(-x)
	null = 0*x
	# plt.figure()
	# plt.fill_between(x, y1, y2)
	plt.plot(x, y, label='b = ' + str(b))
	plt.plot(x, null, "--b")

b_list = [0.2, 0.5, 0.8, 1, 1.2, 1.5, 2] #, 5, 10]
for b in b_list:
	plotter(b)

plt.figure()
for b in b_list:
	stability_plotter(b)
plt.title("Stability")
plt.legend(loc="best")
plt.xlabel('$u(\\tau)$')
plt.ylabel('$du/d\\tau$')
# plt.xscale("log")
plt.savefig("fig/SIR_stability.pdf")

# plot u* as func of b
u_list = []
b_max = 2 
b_list = np.arange(0.01, b_max, 0.01)
for b in b_list:
	func = lambda x : a - b*x - np.exp(-x)
	guess = a/b
	u_eq = fsolve(func, guess)
	u_list.append(u_eq)
plt.figure()
plt.plot(b_list, u_list)
plt.title("equilibrium of u, parametrized by 'b'")
plt.legend(loc="best")
plt.xlabel('$b$')
plt.ylabel('$u*$')
# plt.xscale("log")
plt.savefig("fig/SIR_equilibrium_u_b.pdf")



