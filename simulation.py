import numpy as np

def grav_acc(r, M, G):

    x, y = r
    a_x = - (G*M*x)/(((x**2+y**2)**(1/2))**3)
    a_y = - (G*M*y)/(((x**2+y**2)**(1/2))**3)

    return np.array([a_x, a_y])

def vel_verlet(r_n, v_n, a_n, h ,  M, G):

    r_n1 = r_n + v_n * h + 1/2*a_n * h**2

    a_n1 = grav_acc(r_n1, M, G)

    v_n1 = v_n + 1/2* (a_n + a_n1) * h

    return r_n1, v_n1, a_n1