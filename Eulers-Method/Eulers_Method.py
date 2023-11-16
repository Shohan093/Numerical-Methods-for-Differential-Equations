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
    y = (t + 1) ** 2 - 0.5 * np.exp(t)
    return y

# initial values
t0 = 0.0
y0 = 0.5
t_final = 2
h = 0.1 # step size
n = int((t_final - t0) / h)
t_val = np.zeros(n + 1)
y_val = np.zeros(n + 1)
t_val[0] = t0
y_val[0] = y0

# applying the Euler's Method
for i in range(1, n + 1):
    t_val[i] = round(t_val[i - 1] + h, 3)
    y_val[i] = round(y_val[i - 1] + h * dydt(t_val[i - 1], y_val[i - 1]), 7)

# printing the solution
for i in range(0, n + 1):
    exact = round(exact_sol(t_val[i]), 7)
    print("t: ", t_val[i], "\ty: ", y_val[i], "\t\texact solution: ", exact)

# plotting the solution
plt.figure('Solution')
plt.plot(t_val, y_val, 'r-o')
plt.grid(True)
plt.xlabel('t')
plt.ylabel('y')
plt.legend('N')
plt.title('Numeric Solution to a IVP')
plt.show()
