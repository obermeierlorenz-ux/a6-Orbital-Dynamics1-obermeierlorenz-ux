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

#t for 1 complete circular 

t1 = (2*r0 * math.pi)/ v_c
t10 = 10 * t1


n_steps = 10000
h =  t10 / n_steps

r_0 = np.array([1, 0])


v0_08 = 0.8
v0_12 = 1.2
v_12 = np.array([0, 1.2])
v_08 = np.array([0, 0.8])


a_08 = 1 /(2 /r0 - v0_08**2 /(G * M))
a_12 = 1 /(2 /r0 - v0_12**2 /(G * M))

T_e08 = 2 * math.pi * (a_08**3 /(G *M))**(1/2)
T_e12 = 2 * math.pi * (a_12**3 /(G *M))**(1/2)

print('Expected semimajor axis for v = 0.8:', a_08)
print('Expected period for v = 0.8:', T_e08)
print('Expected semimajor axis for v = 1.2:', a_12)
print('Expected period for v = 1.2:', T_e12)


n_steps = 10000

t_end08 = 5 *T_e08
t_end12 = 5 *T_e12

h08 = t_end08 /n_steps
h12 = t_end12 /n_steps

# simulation for v=0.8

r08 = r_0
v08 = v_08
a08 = grav_acc(r08, M, G)

t_values08 = [0.0]
r_values08 = [r08]
v_values08 = [v08]


for i in range(n_steps):

    r08, v08, a08 = vel_verlet(r08, v08, a08, h08, M, G)

    t_values08.append(t_values08[i] +h08)
    r_values08.append(r08)
    v_values08.append(v08)

t_values08 = np.array(t_values08)
r_values08 = np.array(r_values08)
v_values08 = np.array(v_values08)


# simulation for v01.2

r12 = r_0
v12 = v_12
a12 = grav_acc(r12, M, G)

t_values12 = [0.0]
r_values12 = [r12]
v_values12 = [v12]

for i in range(n_steps):

    r12, v12, a12 = vel_verlet(r12, v12, a12, h12, M, G)

    t_values12.append(t_values12[i] +h12)
    r_values12.append(r12)
    v_values12.append(v12)


t_values12 = np.array(t_values12)
r_values12 = np.array(r_values12)
v_values12 = np.array(v_values12)


# positionen berechnen

x_values08 = r_values08[:, 0]
y_values08 = r_values08[:, 1]
x_values12 = r_values12[:, 0]
y_values12 = r_values12[:, 1]

#Definiton zum berechnen der Energie 
def calculate_energy(r_values,v_values) :

    K_values = []
    U_values = []
    E_values = []

    for i in range(len(r_values)):

        K = 1 / 2 * m * np.sum(v_values[i]**2)

        r_magnitude = np.linalg.norm(r_values[i])

        U = -(G * M * m) / r_magnitude

        E = K + U

        K_values.append(K)
        U_values.append(U)
        E_values.append(E)

    return np.array(K_values), np.array(U_values), np.array(E_values)

#
K_values08, U_values08, E_values08 = calculate_energy(r_values08, v_values08)
K_values12, U_values12, E_values12 = calculate_energy(r_values12, v_values12)


# distance from the central mass and speed

r_magnitude08 = np.linalg.norm(r_values08)
r_magnitude12 = np.linalg.norm(r_values12)

speed08 = np.linalg.norm(v_values08)
speed12 = np.linalg.norm(v_values12)

# measure the period from repeated periapsis locations

def measure_period(t_values,r_magnitude):

    periapsis_times = []

    for i in range(1, len(r_magnitude) - 1):

        if r_magnitude[i] <r_magnitude[i - 1] and r_magnitude[i] <r_magnitude[i + 1]:

            periapsis_times.append(t_values[i])

    periapsis_times = np.array(periapsis_times)

    periods = np.diff(periapsis_times)

    return np.mean(periods)



measured_period08 = measure_period(t_values08, r_magnitude08)
measured_period12 = measure_period(t_values12, r_magnitude12)

fractional_error08 = abs(measured_period08 - T_e08) / T_e08
fractional_error12 = abs(measured_period12 - T_e12) / T_e12


print('v = 0.8')
print('Approximate periapsis distance:', np.min(r_magnitude08))
print('Approximate apoapsis distance:', np.max(r_magnitude08))
print('Measured period:', measured_period08)
print('Fractional period error:', fractional_error08)
print('v = 1.2')
print('Approximate periapsis distance:', np.min(r_magnitude12))
print('Approximate apoapsis distance:', np.max(r_magnitude12))
print('Measured period:', measured_period12)
print('Fractional period error:', fractional_error12)