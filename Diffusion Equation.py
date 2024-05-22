# %%

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters

alpha = 0.01 # Thermal diffusivity
Lx, Ly = 1, 1 # Dimensions of rectangular domain
nx, ny = 50, 50 # NUmber of grid points
dx, dy = (Lx / (nx -1)), (Ly / (ny -1))
T = 10 # Total time (seconds)
dt = 0.01 # Time step

# Creating grid

x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)

X, Y = np.meshgrid(x,y)

# Coefficients for the Crank-Nicolson scheme
rx = alpha * dt / dx**2
ry = alpha * dt / dy**2

a = -0.5 * ry
b = -0.5 * rx
c = 1 + rx + ry


# Boundary conditions
u = np.zeros((nx, ny))

u[:, 0] = 0  # Left boundary
u[:, -1] = 0  # Right boundary
u[0, :] = 0  # Bottom boundary
u[-1, :] = 10  # Top boundary

# Initial conditions

u_initial = np.zeros((nx, ny))

# Creating A matrix

A = np.zeros((nx, ny))

for i in range(nx):
    A[i, i] = c
    if i < nx -1:
        A[i, i+1] = b
    if i > 0:
        A[i, i-1] = b
    if i < (nx - 4): 
        A[i, i + 4] = a
    if i > 3:
        A[i, i - 4] = a

# %%
# Creating d matrix

d = np.zeros((nx, ny))



for i in range(1, nx-1):
    for j in range(1, ny-1):
        d[i][j] = 1 + 0.5 * rx * (u[i+1, j] - 2*u[i, j] + u[i-1, j]) + 0.5 * ry * (u[i, j+1] - 2*u[i, j] + u[i, j-1])


# %%

for i in range(int((T / dt))):
    u = np.linalg.solve(A, d)

print(u)









num_steps = 5
time_points = np.linspace(0, T, num_steps)

for t in time_points:
    # Create d matrix
    d = np.zeros((nx, ny))
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            d[i][j] = 1 + 0.5 * rx * (u[i + 1, j] - 2 * u[i, j] + u[i - 1, j]) + 0.5 * ry * (u[i, j + 1] - 2 * u[i, j] + u[i, j - 1])
    
    # Solve the system of equations
    u = np.linalg.solve(A, d)
    
    # Plot the temperature distribution in 3D
    ax = plt.axes(projection='3d')
    ax.plot3D(X, Y, u)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Temperature')
    plt.title(f'Temperature Distribution at Time {t:.2f} s')
    plt.show()

# # Create figure for 3D plot
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# # Plot the 3D surface
# surf = ax.plot_surface(X, Y, u, cmap='viridis')
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Temperature')
# plt.title('Temperature Distribution')

# # Add a color bar which maps values to colors
# fig.colorbar(surf)

# # Show plot
# plt.show()





# %%
