# -*- coding: utf-8 -*-
"""

Created on Sun September 7 19:38:18 2014

Author: SeanTalia

The code below simulates the currents flowing through the segments of a discrete
electric transmission line.

As of now, the user is free to enter values for the R, L, G, and C in each segment
as well as the number of segments into which the line has been divided, although
the code will only plot up to the first 35 segments. The user is also free to
play around with the resistance and inductance of the load.

The code works as follows. We have a voltage function V(t) = Acos(omega*t) which
applies a voltage to a line made up of N segments. Segment j of the line has a
current I_j flowing through it as well as a voltage V_j associated with it. If
there are N segments, we can use Ohm's Law on each segment to end up with an
algebraic system of 2N+1 equations:

    V - V_1 = ZI_1
    V_1 - V_2 = ZI_2
    .
    .
    .
    V_{n-1} - V_n = ZI_n
    ------------------------
    I_1 - I_2 = YV_1
    I_2 - I_3 = YV_2
    .
    .
    .
    I_n - I_load = YV_n
    ------------------------
    V_n = Z_load*I_load

The V_j's and the I_j's are the complex coefficients on expressions of the form
I_j*exp(i*omega*t), so the complex coefficient I_j determines the amplitude and
phase of the current passing through segment j.

The function below determines the complex coefficient for each V_j and I_j by
solving this algebraic system. The user may choose to plot the currents up to
the first 35 stations (or more, if they choose to add more colors to the
colors array on line 159.)

"""
import numpy as np
import matplotlib.pyplot as plt

# Current accepts a complex number Z which determines the amplitude and phase
# and returns the desired sinusoidal function. This takes a complex-valued
# sinusoidal function Ze^(i*omega*t), and returns the real-valued sinusoidal
# function that corresponds to the real part of the complex-valued function.

def current(Z, t, omega):

    # Amplitude
    A = np.abs(Z)

    x = Z.real
    y = Z.imag

    # Phase
    theta = np.arctan2(y, x)

    return A*np.cos(omega*t + theta)

# The function below does all the heavy lifting - it accepts the series resistance
# and inductance, the shunt admittance and capacitance, the load resistance and
# inductance, the number of transmission line segments, and the voltage
# amplitude and frequency as parameters.

def transLine(R, L, G, C, R_l, L_l, m, omega, V):

    # System Dimension
    N = 2*m + 1

    # Series Impedance and Shunt Admittance Per Segment
    Z = (R + omega*L*(1.0j)) / m
    Y = (G + omega*C*(1.0j)) / m

    # Load Impedance
    Z_l = R_l + omega*L_l*(1.0j)

    # Construct System Using 3 Parts
    voltages = np.zeros((m, N), dtype = complex)
    currents = np.zeros((m, N), dtype = complex)
    load_eq  = np.zeros((N,), dtype = complex)

    # Modeling Voltage Equations
    for i in range(m):
      # 1's on diagonal
      voltages[i][i] = 1.0
      # -1's to the left of the 1's
      if i > 0:
          voltages[i][i-1] = -1.0
      # Z's on (a different) diagonal
      voltages[i][m+i] = Z

    # Modeling Current Equations
    for i in range(m):
      # Y's on diagonal
      currents[i][i] = Y
      # -1's on (a different) diagonal
      currents[i][m+i] = -1.0
      # 1's to the right of -1's
      currents[i][m+i+1] = 1.0

    # Modeling Load Equation (V_n = Z_l*I_l)
    load_eq[m-1] = -1
    load_eq[N-1] = Z_l

    # Constructing our matrix of V's and I's
    A = np.vstack([np.concatenate((voltages, currents)), load_eq])

    # Solution Vector
    b = np.zeros((N, ), dtype = complex)
    b[0] = V
    b[N-1] = 0

    # Solving the System
    x = np.linalg.solve(A, b)

    # Checking that the solution is correct
    print np.allclose(np.dot(A, x), b)

    # Collecting the coefficients for our current functions
    I = np.zeros((m, ), dtype = complex)

    for i in range(m):
      I[i] = x[m+i]

    # Return vector of complex coefficients
    return I

# This function accepts an array (of up to length 35) I of complex coefficients and plots the
# sinusoidal curves associated with each of these coefficients.
def currentPlot(I, T, omega, m):

    # Number of points to plot
    p = 5001

    # Equispaced points for plotting
    x_vals = np.zeros(p)

    for i in range(p):
      x_vals[i] = i * (T / (p-1))

    # The remainder of the function creates a plot and graphs the desired number of curves
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # Labeling Axes
    ax.set_xlabel('Time')
    ax.set_ylabel('Current')

    # Curve Colors
    colors = ['#ff1493', '#ff0000', '#ff4500', '#ff6347', '#ff7f50',
              '#ffa500', '#ffb600', '#ffd700', '#fff000', '#ffff00',
              '#9acd32', '#adff2f', '#00ff00', '#00cc44', '#00aa22',
              '#228b22', '#7fff7f', '#4cff4c', '#00ffff', '#00fa9a',
              '#00bfff', '#4169e1', '#0000cd', '#000080', '#191970',
              '#5518ab', '#4c177d', '#44146f', '#3b1261', '#330f53',
              '#2a0d45', '#220a37', '#190729', '#11051b', '#82020d']

    # Plotting our curves
    for i in range(m):
      ax.plot(x_vals, current(I[i], x_vals, omega), color = colors[i])

    return plt

def main():

    # Voltage Frequency and Amplitude
    omega = 0.4
    V = 5.0

    # Total series resistance and inductance, shunt admittance and capacitance of whole line
    R = 0.8
    L = 3.5
    G = 0.0
    C = 2.5

    # Load Resistance and Inductance
    R_l = 10.0
    L_l = 15.0

    # Number of segments into which line has been divided
    N = 15

    # Final Time
    T = 80.0

    # Get our array of complex coefficients for the currents
    I = transLine(R, L, G, C, R_l, L_l, N, omega, V)

    # Plot the currents
    currentPlot(I, T, omega, N).show()

main()
