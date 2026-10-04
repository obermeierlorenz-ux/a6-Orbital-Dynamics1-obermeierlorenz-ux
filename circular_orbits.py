import numpy as np
import math
from simulation import grav_acc, vel_verlet
import matplotlib.pyplot as plt

G = 1
M = 1
m = 1
r0 = 1
v_c= 1
T_c = 2*math.pi

#t for 1 compete circular 

t1 = (2*r0 * math.pi)/ v_c
t10 = 10 * t1


n_steps = 10000
h =  t10 / n_steps

r_0 = np.array([1, 0])
v_0 = np.array([0, 1])

a_0= grav_acc(r_0, M, G)

t_values = [0.0]
r_values = [r_0]
v_values = [v_0]
a_values = [a_0]

r = r_0
v = v_0
a = a_0

for i in range (n_steps):#

    r, v, a = vel_verlet(r, v, a, h, M, G)

    t_values.append(t_values[i] + h)
    r_values.append(r)
    v_values.append(v)
    a_values.append(a)



t_values = np.array(t_values)
r_values = np.array(r_values)
v_values = np.array(v_values)
a_values = np.array(a_values)

x_values = r_values[:, 0]
y_values = r_values[:, 1]

fig1 = plt.figure()
ax1 = fig1.add_axes([0.12, 0.12, 0.80, 0.80])
ax1.plot(0, 0, marker='*', color='orange',
         markersize=12, label='Central mass')
ax1.plot(x_values, y_values, label = 'Orbits')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('Circular Orbit Trajectory')
# same visible scale for x and y
ax1.set_aspect('equal', adjustable='box')
ax1.legend()
ax1.grid()
plt.show()

fig2 = plt.figure()
ax2 = fig2.add_axes([0.12, 0.12, 0.80, 0.80])
ax2.plot(t_values, x_values, label = 'x(t)')
ax2.plot(t_values, y_values, label = 'y(t)')
ax2.set_xlim(0, t10)
ax2.set_xlabel('Time t')
ax2.set_ylabel('Position')
ax2.set_title('Position Components over Ten Circular Orbits')
ax2.grid()
ax2.legend()
plt.show()

#Energie vs time 

K_values = []
U_values = []
E_values = []

for iters in range(len(t_values)): 


    K = 1/2 * M * np.sum(v_values[iters]**2)

    # gibt mir den betreag des r Vektors
    r_magnitude = np.linalg.norm(r_values[iters])

    U = - (G*M*m) / r_magnitude
    E = K + U 

    K_values.append(K)
    U_values.append(U)
    E_values.append(E)

K_values = np.array(K_values)
U_values = np.array(U_values)
E_values = np.array(E_values)


fig3 = plt.figure()
ax3 = fig3.add_axes([0.12, 0.12, 0.80, 0.80])
ax3.plot(t_values, K_values, label = 'kinetic energy')
ax3.plot(t_values, U_values, label = 'potential energy')
ax3.plot(t_values, E_values, label = 'total mechanical ernergy')
ax3.set_xlabel('Time t')
ax3.set_ylabel('Energy')
ax3.set_title('K(t), U(t) and E(t)')
ax3.grid()
ax3.legend()
plt.show()

#period from trajectory bei looking at the y values

crossing_times = []

for i in range(1, len(y_values)):

    if y_values[i-1] < 0 and y_values [i] >=0: 
        crossing_time = t_values[i]
        crossing_times.append(crossing_time)

periods = []

for i in range(len(crossing_times)-1): 

    period = crossing_times[i+1] - crossing_times[i]

    periods.append(period)

#print(periods)
periods = np.array(periods)
mean_periods = np.mean(periods)

periods_diff = T_c - mean_periods

print(periods_diff)
