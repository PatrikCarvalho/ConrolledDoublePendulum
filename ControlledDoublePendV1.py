"""
===========================
The double pendulum problem
===========================

This animation illustrates the double pendulum problem.

Double pendulum formula translated from the C code at
http://www.physics.usyd.edu.au/~wheat/dpend_html/solve_dpend.c
"""

from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
import A_Lin
import control.matlab as ctr
from matplotlib.gridspec import GridSpec

#global K, state_eql, Control

G = 9.81  # acceleration due to gravity, in m/s^2
L1 = 1 # length of pendulum 1 in m
L2 = 1  # length of pendulum 2 in m
M1 = 2   # mass of pendulum 1 in kg
M2 = 1   # mass of pendulum 2 in kg
T_1 = 4  # torque on the first hinge [Nm]

# Linearized A and B matrices, must not confuse equilibrium point with initial point
th1_eql = 180
w1_eql = 0
th2_eql = 180
w2_eql = 0
state_eql = np.radians([th1_eql, w1_eql, th2_eql, w2_eql])

# Ofset from the Equilibrium point in degrees for the initial conditions
th1_offset_eql = 30
w1_offset_eql = 0
th2_offset_eql = -30
w2_offset_eql = 0

A_lin = np.matrix(A_Lin.A_lin(M1, M2, L1, L2, G, state_eql[0], state_eql[1], state_eql[2], state_eql[3]), dtype=np.float)
B_lin = np.matrix([0, 1, 0, 0]).T
C_lin = np.matrix(np.identity(len(A_lin)))
D_lin = np.matrix(np.zeros_like(B_lin))

#sys = ctr.ss(A_lin, B_lin, C_lin, D_lin)
sys = ctr.StateSpace(A_lin, B_lin, C_lin, D_lin)
Q = np.matrix([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1000, 0],
               [0, 0, 0, 1000]])
R = 1
K, S, E = ctr.lqr(sys, Q, R)




def derivs(state, t):
    if Control == True:
        T_1 = -K.dot(state - state_eql)
        T_1 = T_1[0]
    elif Control == False:
        T_1 = 0

    dydx = np.zeros_like(state)
    dydx[0] = state[1]
    
    delta = state[2] - state[0]
    den1 = (M1+M2) * L1 - M2 * L1 * cos(delta) * cos(delta)
    dydx[1] = ((M2 * L1 * state[1] * state[1] * sin(delta) * cos(delta) + M2 * G * sin(state[2]) * cos(delta) + M2 * L2 * state[3] * state[3] * sin(delta) - (M1+M2) * G * sin(state[0]) + T_1) / den1)

    dydx[2] = state[3]
    
    den2 = (L2/L1) * den1

    dydx[3] = ((- M2 * L2 * state[3] * state[3] * sin(delta) * cos(delta) + (M1+M2) * G * sin(state[0]) * cos(delta) - (M1+M2) * L1 * state[1] * state[1] * sin(delta) - (M1+M2) * G * sin(state[2])) / den2)

    return dydx

# create a time array from 0..100 sampled at 0.05 second steps
dt = 0.05
t = np.arange(0, 10, dt)


# th1 and th2 are the initial angles (degrees)
# w10 and w20 are the initial angular velocities (degrees per second)
th1 = th1_eql + th1_offset_eql
w1 = w1_eql + w1_offset_eql
th2 = th2_eql + th2_offset_eql
w2 = w2_eql + w2_offset_eql

# initial state
state_control = np.radians([th1, w1, th2, w2])
state_not_control = state_control

# integrate your ODE using scipy.integrate.
Control = True
y_control = integrate.odeint(derivs, state_control, t)
T1_control = np.squeeze(-K.dot(y_control.transpose()))

Control = False
y_not_control = integrate.odeint(derivs, state_not_control, t)
T1_not_control = np.zeros(np.shape(y_not_control[:,0]))

x1_control = L1*sin(y_control[:, 0])
y1_control = -L1*cos(y_control[:, 0])
x1_not_control = L1*sin(y_not_control[:, 0])
y1_not_control = -L1*cos(y_not_control[:, 0])

x2_control = L2*sin(y_control[:, 2]) + x1_control
y2_control = -L2*cos(y_control[:, 2]) + y1_control
x2_not_control = L2*sin(y_not_control[:, 2]) + x1_not_control
y2_not_control = -L2*cos(y_not_control[:, 2]) + y1_not_control

fig = plt.figure()
gs = GridSpec(2, 2, figure=fig)
ax1 = fig.add_subplot(221, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
ax1.set_aspect('equal')
ax1.grid()
ax1.set_title("Double Pendulum with Control")
ax2 = fig.add_subplot(222, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
ax2.set_aspect('equal')
ax2.grid()
ax2.set_title("Double Pendulum without Control")
ax3 = fig.add_subplot(gs[1, :], autoscale_on=False, xlim=(0, np.max(t)), ylim=(np.min(T1_control), np.max(T1_control)))
ax3.grid()
ax3.set_title("Applied Torque")

line1, = ax1.plot([], [], 'o-', lw=2)
line2, = ax2.plot([], [], 'o-', lw=2)
line3, = ax3.plot([], [], lw=2)
time_template = 'time = %.1fs'
time_text1 = ax1.text(0.05, 0.9, '', transform=ax1.transAxes)
time_text2 = ax2.text(0.05, 0.9, '', transform=ax2.transAxes)
time_text3 = ax3.text(0.05, 0.9, '', transform=ax3.transAxes)


def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    time_text1.set_text('')
    time_text2.set_text('')
    time_text3.set_text('')
    return line1, time_text1, line2, time_text2, line3, time_text3


def animate(i):
    thisx_control = [0, x1_control[i], x2_control[i]]
    thisy_control = [0, y1_control[i], y2_control[i]]
    thisx_not_control = [0, x1_not_control[i], x2_not_control[i]]
    thisy_not_control = [0, y1_not_control[i], y2_not_control[i]]
    thist = t[0:i+1]
    thisT1_c = T1_control[0:i+1]


    line1.set_data(thisx_control, thisy_control)
    time_text1.set_text(time_template % (i*dt))
    line2.set_data(thisx_not_control, thisy_not_control)
    time_text2.set_text(time_template % (i*dt))
    line3.set_data(thist, thisT1_c)
    time_text3.set_text(time_template % (i*dt))
    return line1, time_text1, line2, time_text2, line3, time_text3


ani = animation.FuncAnimation(fig, animate, range(1, len(y_control)),
                              interval=dt*1000, blit=True, init_func=init)
plt.show()
