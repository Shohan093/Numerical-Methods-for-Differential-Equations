import numpy as np
import matplotlib.pyplot as plt

# Parameters
Lx, Ly = 10.0, 10.0  # Domain size
Nx, Ny = 50, 50      # Number of grid points
dx, dy = Lx/(Nx-1), Ly/(Ny-1)  # Grid spacing

# Initial guess for the solution
u = np.zeros((Ny, Nx))

# Boundary conditions
u[:, 0] = 100  # Left boundary
u[:, -1] = 100  # Right boundary
u[0, :] = 100  # Bottom boundary
u[-1, :] = 100  # Top boundary

# Iterative solver parameters
tolerance = 1e-5
max_iterations = 10000

# Solving the Laplace equation using the finite difference method
for iteration in range(max_iterations):
    u_old = u.copy()
    
    # Update interior points
    u[1:-1, 1:-1] = 0.25 * (u_old[1:-1, :-2] + u_old[1:-1, 2:] + u_old[:-2, 1:-1] + u_old[2:, 1:-1])
    
    # Check for convergence
    if np.linalg.norm(u - u_old) < tolerance:
        break

# Plotting the solution
X, Y = np.meshgrid(np.linspace(0, Lx, Nx), np.linspace(0, Ly, Ny))
plt.contourf(X, Y, u, cmap='hot')
plt.colorbar(label='Temperature')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Solution of the Laplace Equation')
plt.show()
