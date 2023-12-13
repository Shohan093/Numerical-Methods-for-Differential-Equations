# solution to the differential equation y' = y - t^2 + 1, y(0) = 0.5, with 0<=t<=2
# exact solution for this differential equaiton is y = (t + 1)^2 - 0.5 * exp(t)
# for more accurate results decrease the step size

from typing import Any
import numpy as np
from matplotlib import pyplot as plt

# defining the IVP
def dydt(t, y):
    y1: Any = y - t ** 2 + 1
    return y1

# exact solution
def exact_sol(t):
    y: Any = (t + 1) ** 2 - 0.5 * np.exp(t)
    return y

# Runge-Kutta of second order
def rk2(f, t0, t_final, y0, h):
    n = int((t_final - t0) / h)
    t_val = np.zeros(n + 1)
    y_val = np.zeros(n + 1)
    t_val[0] = t0
    y_val[0] = y0
    
    # applying the second order Runge Kutta Method
    for i in range(1, n + 1):
        t_val[i] = round(t_val[i - 1] + h, 3)
        k1 = h * f(t_val[i - 1], y_val[i - 1])
        k2 = h * f(t_val[i - 1] + h, y_val[i - 1] + k1)
        y_val[i] = round(y_val[i - 1] + 0.5 * (k1 + k2), 7)
    return t_val, y_val
    

# initial values
t0 = 0.0
y0 = 0.5
t_final = 2
h = 0.1 # step size
t, y = rk2(dydt, t0, t_final, y0, h)

# printing the iterations
n = int((t_final - t0) / h)
for i in range(0, n + 1):
    exact = round(exact_sol(t[i]), 7)
    print('t: ', t[i], '\ty: ', y[i], '\texact solution: ', exact)

# ploting the solution
plt.figure('Solution', figsize=(10, 8))
plt.plot(t, y, 'b-o')
plt.xlabel('t')
plt.ylabel('y')
plt.grid(True)
plt.title('Numerical solution using RK2')
plt.show()
