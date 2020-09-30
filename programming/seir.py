import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import os


data = {}
directory = "contact_matrices"
for filename in os.listdir(directory):
    if filename.endswith(".csv"): 
        # print(os.path.join(directory, filename))
        fn = filename.split(".")[0]
        filename_path = os.path.join(directory, filename)
        data[fn] = np.genfromtxt(filename_path, delimiter=',')
        # Uncomment the following two lines for matrix visualizasion.
        plt.imshow(data[fn], cmap='hot', interpolation='nearest')
        # plt.show()
        continue

N = 1e6
i0 = 1 #initial infection
R0 = 2.68 #basic reproduction number
a = 1/3. #reciprocal incubation time.
# b = 1/3. #reciprocal time between contacts
g = 1/4. #reciprocal recovery time 
eta = 0.26 # infection efficiency

# function that returns dP/dt
def model(P,t):
    S = P[0] # Susceptible
    E = P[1] # Exposed
    I = P[2] # Infected
    R = P[3] # removed
    dSdt = - R0*g*I*S/N 
    dEdt = R0*g*I*S/N - a*E 
    dIdt = a*E - g*I 
    dRdt = g*I
    dPdt = [dSdt,dEdt,dIdt,dRdt]
    return dPdt

# initial condition
P0 = [N - i0, 0, i0, 0]

# time points
t = np.linspace(0,int(10))

# solve ODE
P = odeint(model,P0,t)

# plot results
plt.plot(t,P)
plt.legend(["Susceptible", "Exposed", "Infected", "Removed" ], loc="best")
plt.xlabel('time')
plt.ylabel('P(t)')
plt.savefig("test.pdf")
# plt.show()
