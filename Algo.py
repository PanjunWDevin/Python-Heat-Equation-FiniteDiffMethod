# Heat Equation: Explicit Method

from matplotlib import pylab
import numpy as np
import math as mt
from matplotlib import pylab

#build initial conditions
def initial_con(x):
    return 4 * x * (1 - x )

#print initial condiiton image
#nodes_x_ini = np.linspace(0,1,50)
#nodes_u_ini = map(initial_con,nodes_x_ini)

#pylab.plot(nodes_x_ini,nodes_u_ini)
#pylab.show()

#Set Parameter for calculating Heat Equations u_t = k * u_xx by forward finite difference method

N = 25    # x axis number of steps
M = 2500  # t axis number of steps

X = 1.0   # X axis Bound
T = 1.0   # T axis Bound

u_nodes = np.zeros([M+1,N+1]) # boundary conditions set at x = 0, nodes = 0, at x = 1, nodes = 0
u_nodes[0,:] = np.linspace(0,1,N+1)
xArray = np.linspace(0,1,N+1)
u_nodes[0,:] = map(initial_con,u_nodes[0,:]) #assign the initial conditions as u(x,0) = 4 * x * (1-x)
#print u_nodes[0,:]

x_step = X/N #initialise the x step
t_step = T/M #initialise the time step

#calculate the node value in the matrix node by finite difference method
#set parameter pho

kappa = 1.0 #coffience for head transformation

pho = kappa * t_step / (x_step * x_step)
#print pho

for k in range(1,M+1,1): #only run till M, namely stops at index = M-1, cause we have to stop at t = T
    for j in range (1,N,1):
        u_nodes[k][j] = pho * u_nodes[k-1][j-1] + (1 - 2 * pho) * u_nodes[k-1][j] + pho * u_nodes[k-1][j+1]

#print u_nodes[M,:] #test output

#Plot
pylab.figure(figsize = (12,6))
pylab.plot(xArray, u_nodes[0,:])
pylab.plot(xArray, u_nodes[int(0.10/ t_step),:])
pylab.plot(xArray, u_nodes[int(0.20/ t_step),:])
pylab.plot(xArray, u_nodes[int(0.50/ t_step),:])
pylab.xlabel('$x$', fontsize = 15)
pylab.ylabel(r'$U(\dot, \tau)$', fontsize = 15)
pylab.title('one dimensional heat equation')
pylab.legend([r'$\tau = 0.$', r'$\tau = 0.10$', r'$\tau = 0.20$', r'$\tau = 0.50$'], fontsize = 15)
pylab.show()

#3-D Plot about the results
