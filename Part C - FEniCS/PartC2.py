# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 23:44:48 2021

@author: matil

#This code is used in Part C2
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import loadtxt
from dolfin import *
 
# Create mesh and define function space
mesh = Mesh("circle_fine.xml.gz")

# Construct the finite element space
V = VectorFunctionSpace (mesh, 'P', 1)

# Define parameters
T = 1000
k = 0.5
alpha = 0.4
beta = 2
gamma = 0.8
delta1 = 1
delta2 = 1

   
# Class representing the intial conditions
# Randomly generated values for problem C2
class InitialConditions(UserExpression):
    def eval(self, values,  x):
        values[0] = 1/2*(1-np.random.uniform(0,1))
        values[1] = 1/4+1/2*np.random.uniform(0,1)
        
    def value_shape(self):
        return (2,)


# Define initial condition
indata = InitialConditions(degree=2)
u_ini = Function(V)
u_ini = interpolate(indata,V)


# Test and trial functions
u = TrialFunction(V)
v = TestFunction(V)

# Create bilinear and linear forms
a0 = u[0]*v[0]*dx -1/2*k*u[0]*v[0]*dx + 1/2*k*delta1*inner(grad(u[0]), grad(v[0]))*dx  
a1 = u[1]*v[1]*dx + gamma*k*1/2*u[1]*v[1]*dx + 1/2*k*delta2*inner(grad(u[1]), grad(v[1]))*dx 
      
L0 = u_ini[0]*v[0]*dx + 1/2*k*u_ini[0]*v[0]*dx -\
    (1/2*k*delta1*inner(grad(u_ini[0]), grad(v[0]))*dx) -\
    (k*((u_ini[0]*u_ini[1])/(u_ini[0]+alpha) + u_ini[0]*u_ini[0])*v[0]*dx)  
    
L1 = u_ini[1]*v[1]*dx - gamma*1/2*k*u_ini [1]*v[1]*dx -\
    (1/2*k*delta2*inner(grad(u_ini[1]), grad(v[1]))*dx) -\
    (k*(-(beta*u_ini[1]*u_ini[0])/(u_ini[0]+alpha))*v[1]*dx)  
    
a = a0 + a1
L = L0 + L1

# Set an output file
file = File("C2_fine/sol_c2.pvd","compressed")

# Calculating solution
u = Function(V)
u.assign(u_ini)
pop_prey = []
pop_pred = []
time = []
M_prey = u[0]*dx
M_pred = u[1]*dx
prey = assemble(M_prey)
pred = assemble(M_pred)

# Time step
t = 0.0
file << (u,t)
while t <= T:
    time.append(t)
    u_ini.assign(u)
    A = assemble(a)
    b = assemble(L)
    
    solve (A, u.vector(), b, "lu")
    
    prey = assemble(M_prey)
    pred = assemble(M_pred)
    pop_prey.append(prey)
    pop_pred.append(pred)
    
    #Saving solution to file in steps of t = 50
    if t % 50 == 0:
        print("Saved, t = " + str(t))
        file << (u,t)
        
    t += k
    

# Plotting stuff
plt.plot(time,pop_prey, label = 'Prey')
plt.plot(time, pop_pred, label = 'Predators')
plt.title("Prey/Predator Dynamics")
plt.xlabel("Time, t")
plt.ylabel("Amount, #")
plt.legend(loc = 'upper right')

# Saving stuff
plt.savefig("C2_fine/plot.png")
np.savetxt("C2_fine/time.txt", time)
np.savetxt("C2_fine/pop_prey.txt", pop_prey)
np.savetxt("C2_fine/pop_pred.txt", pop_pred)

