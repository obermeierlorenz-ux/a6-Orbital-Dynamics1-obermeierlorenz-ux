import numpy as np
import math
from simulation import grav_acc, vel_verlet
import matplotlib.pyplot as plt

G = 1
M = 1
m = 1
r0 = 1
v_c = (G * M / r0)**(1 / 2)
v_escape = (2 * G * M / r0)**(1 / 2)

T_c = 2*math.pi

#t for 1 compete circular 

t1 = (2*r0 * math.pi)/ v_c
t10 = 10 * t1


n_steps = 10000
h =  t10 / n_steps

r_0 = np.array([1, 0])

v_12 = np.array([0, 1.2])
v_15 = np.array([0, 1.5])

a_0= grav_acc(r_0, M, G)

t_values = [0.0]
r_values = [r_0]
v_values12 = [v_12]
a_values = [a_0]

r = r_0
v12 = v_12
a = a_0

#v12
for i in range (n_steps):#

    r, v12, a = vel_verlet(r, v12, a, h, M, G)

    t_values.append(t_values[i] + h)
    r_values.append(r)
    v_values12.append(v12)
    a_values.append(a)


t_values12 = np.array(t_values)
r_values12 = np.array(r_values)
v_values12 = np.array(v_values12)
a_values12 = np.array(a_values)

#v15

r15 = r_0
v15 = v_15
a15 = grav_acc(r15, M, G)

r_values15 = [r15]
v_values15 = [v15]

for i in range(n_steps):

    r15, v15, a15 = vel_verlet(r15, v15, a15, h, M, G)

    r_values15.append(r15)
    v_values15.append(v15)

r_values15 = np.array(r_values15)
v_values15 = np.array(v_values15)

x_values12 = r_values12[:, 0]
y_values12 = r_values12[:, 1]

x_values15 = r_values15[:, 0]
y_values15 = r_values15[:, 1]


#Referenkreis
theta = np.linspace(0, 2 * math.pi, 300)

circle_x = np.cos(theta)
circle_y = np.sin(theta)

#Plot

fig1 = plt.figure()
ax1 = fig1.add_axes([0.12, 0.12, 0.80, 0.80])
ax1.set_xlim(-3, 3)
ax1.set_ylim(-3, 3)


ax1.plot(0, 0, marker='*', color='orange',
         markersize=12, label='Central mass')
ax1.plot(x_values12, y_values12, label = 'Inital speed v = 1.2')
ax1.plot(x_values15, y_values15, label = 'Initial speed v = 1.5')
ax1.plot(circle_x, circle_y, linestyle='--',
         color='black', label='Reference circle: r = 1')

ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('Close-up Comparison of Orbital Trajectories')

# same visible scale for x and y
ax1.set_aspect('equal', adjustable='box')
ax1.legend()
ax1.grid()
plt.show()


#Energie vs time 

K_values12 = []
U_values12 = []
E_values12 = []

for iters in range(len(t_values12)): 


    K = 1/2 * m * np.sum(v_values12[iters]**2)

    # gibt mir den betreag des r Vektors
    r_magnitude = np.linalg.norm(r_values12[iters])

    U = - (G*M*m) / r_magnitude
    E = K + U 

    K_values12.append(K)
    U_values12.append(U)
    E_values12.append(E)

K_values12 = np.array(K_values12)
U_values12 = np.array(U_values12)
E_values12 = np.array(E_values12)

#energie plot 15
K_values15 = []
U_values15 = []
E_values15 = []

for iters in range(len(t_values)): 


    K = 1/2 * m * np.sum(v_values15[iters]**2)

    # gibt mir den betreag des r Vektors
    r_magnitude = np.linalg.norm(r_values15[iters])

    U = - (G*M*m) / r_magnitude
    E = K + U 

    K_values15.append(K)
    U_values15.append(U)
    E_values15.append(E)

K_values15 = np.array(K_values15)
U_values15 = np.array(U_values15)
E_values15 = np.array(E_values15)

#plot für 1.2 energie
fig2 = plt.figure()
ax2 = fig2.add_axes([0.12, 0.12, 0.80, 0.80])

ax2.plot(t_values12, K_values12, label = 'kinetic energy K')
ax2.plot(t_values12, U_values12, label = 'potential energy U')
ax2.plot(t_values12, E_values12, label = 'total mechanical ernergy E')
ax2.set_xlabel('Time t')
ax2.set_ylabel('Energy')
ax2.set_title('Energy for Initial Speed v = 1.2')
ax2.grid()
ax2.legend()

plt.show()

#plot für 1.5 erngerie
fig3 = plt.figure()
ax3 = fig3.add_axes([0.12, 0.12, 0.80, 0.80])

ax3.plot(t_values12, K_values15, label='Kinetic energy K')
ax3.plot(t_values12, U_values15, label='Potential energy U')
ax3.plot(t_values12, E_values15, label='Total energy E')
ax3.set_xlabel('Time t')
ax3.set_ylabel('Energy')
ax3.set_title('Energy for Initial Speed v = 1.5')
ax3.grid()
ax3.legend()

plt.show()