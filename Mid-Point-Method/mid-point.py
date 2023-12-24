# solution to the differential equation y' = y - t^2 + 1, y(0) = 0.5, with 0<=t<=2
# exact solution for this differential equaiton is y = (t + 1)^2 - 0.5 * exp(t)
# for more accurate results decrease the step size

from typing import Any
from matplotlib import pyplot as plt
import numpy as np

# defining the IVP
def dydt(t, y):
    y1: Any = y - t ** 2 + 1
    return y1

def exact_solution(t):
    y: Any= (t + 1) ** 2 - 0.5 * np.exp(t)
    return y

def mid_point(f, t0, t_final, y0, h):
    n = int((t_final - t0) / h)
    t_val = np.zeros(n + 1)
    y_val = np.zeros(n + 1)
    t_val[0] = t0
    y_val[0] = y0

    # applying the mid point method
    for i in range(1, n + 1):
        t_val[i] = round(t_val[i - 1] + h, 3)
        slope = f(t_val[i - 1], y_val[i - 1]) # next solpe
        est_yval = y_val[i - 1] + 0.5 * h * slope # estimated next y_val
        y_val[i] = round(y_val[i - 1] + h * f(t_val[i - 1] + h * 0.5, est_yval), 7) # next y_val
    return t_val, y_val

if __name__ == "__main__":
    # initial values
    y0 = 0.5
    t0 = 0.0
    t_final = 2
    h = 0.1 # step size
    t, y = mid_point(dydt, t0, t_final, y0, h)

    # printing the estimations
    for i in range(int((t_final - t0) / h) + 1):
        exact = round(exact_solution(t[i]), 7)
        print('t: ', t[i], '\ty: ', y[i], '\texact solution: ', exact)

    # plotting the solutions
    plt.figure(figsize=(8, 6))
    plt.plot(t, y, 'b-o')
    plt.xlabel('t')
    plt.ylabel('y')
    plt.grid(True)
    plt.title('Numerical approximation with Mid-Point Method')
    plt.show()
