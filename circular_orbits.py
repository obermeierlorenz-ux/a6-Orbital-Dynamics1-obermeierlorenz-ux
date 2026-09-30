import numpy as np
import math
from simulation import grav_acc, vel_verlet
import matplotlib.pyplot as plt

G = 1
M = 1
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

    r, v, a = vel_verlet(r, v, a)

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

ax1.plt(x_values, y_values, label = 'Orbits')
plt.plot()
plt.show()