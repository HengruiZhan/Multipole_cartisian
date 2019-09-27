import analytical_potential as ap
import grid
import matplotlib.pyplot as plt
import multipole
import numpy as np
import L2_difference as L2
# from gekko import GEKKO

"""
nx = 64
ny = 64
nz = 64
"""
nx = 128
ny = 128
nz = 128


xlim = (-1.0, 1.0)
ylim = (-1.0, 1.0)
zlim = (-1.0, 1.0)
g = grid.Grid(nx, ny, nz, xlim, ylim, zlim)

dens = g.scratch_array()

# density of the cube
a = 0.3
b = 0.3
c = 0.3
sites = (a, b, c)
# density of a cube
for i in range(g.nx):
    for j in range(g.ny):
        for k in range(g.nz):
            if g.x[i] >= 0 and \
                    g.x[i] <= a and \
                    g.y[j] >= 0 and \
                    g.y[j] <= b and \
                    g.z[k] >= 0 and \
                    g.z[k] <= c:
                dens[i, j, k] = 1.0

# analytical potential of a cube
V_cube = g.scratch_y_plane_array()
for i in range(g.nx):
    for j in range(g.nz):
        cube_V = ap.Analytical_Cube_Potential(sites, g.x[i], g.y[0], g.z[j])
        V_cube[i, j] = cube_V.V
"""
plt.imshow(np.log10(np.abs(np.transpose(V_cube))), origin="lower",
           interpolation="nearest",
           extent=[g.xlim[0], g.xlim[1],
                   g.zlim[0], g.zlim[1]])

plt.colorbar()
ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("V_cube.png")
"""

"""
# plot the analytical potential of the cube
plt.imshow(np.log10(np.abs(np.transpose(V_cube))), origin="lower",
           interpolation="nearest",
           extent=[g.xlim[0], g.xlim[1],
                   g.zlim[0], g.zlim[1]])

plt.colorbar()
ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("phi.png")
"""
"""
# density of a MacLaurin spheroid

sph_center = (0.5, 0.5, 0.5)
a_1 = 0.23
a_3 = 0.10
mask = (g.x3d-sph_center[0])**2/a_1**2 +\
    (g.y3d-sph_center[1])**2/a_1**2 +\
    (g.z3d-sph_center[2])**2/a_3**2 <= 1
dens[mask] = 1.0
density = 1.0
"""


# analytical potential of MacLaurin spheroid
"""
n_moments = 0
# center = (a/2, b/2, c/2)
center = (0.5, 0.5, 0.5)
# multipole expansion with concentric sphere
moment = multipole.Multipole(g, dens, n_moments, 2*g.dr, center=center)
'''
# with azumuthal
for l in range(n_moments):
    for m in range(-l, l+1):
        moment.compute_expansion(dens, l, m)
        '''

phi = g.y_plane_array()

for i in range(g.nx):
    for j in range(g.nz):
        phi[i, j] = moment.phi(g.x[i], g.y[0], g.z[j])
        # phi[i, j] = moment.phi_point(g.x[i], g.y[0], g.z[j])

plt.imshow(np.log10(np.abs(np.transpose(phi))), origin="lower",
           interpolation="nearest",
           extent=[g.xlim[0], g.xlim[1],
                   g.zlim[0], g.zlim[1]])

plt.colorbar()
ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("lmax=0.png")

L2normerr = L2.L2_diff(V_cube, phi, g.scratch_array())
print("resolution", nx)
print("lmax =", n_moments)
print(L2normerr)
"""
# test brute force codes
n_moments = 0
# center = (a/2, b/2, c/2)
center = (0.5, 0.5, 0.5)
# multipole expansion with concentric sphere
moment = multipole.Multipole(g, dens, n_moments, 2*g.dr, center=center)
'''
# with azumuthal
for l in range(n_moments):
    for m in range(-l, l+1):
        moment.compute_expansion(dens, l, m)
        '''

phi = g.y_plane_array()

for i in range(g.nx):
    for j in range(g.nz):
        phi[i, j] = moment.phi(g.x[i], g.y[0], g.z[j])
        # phi[i, j] = moment.phi_point(g.x[i], g.y[0], g.z[j])

plt.imshow(np.log10(np.abs(np.transpose(phi))), origin="lower",
           interpolation="nearest",
           extent=[g.xlim[0], g.xlim[1],
                   g.zlim[0], g.zlim[1]])

plt.colorbar()
ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("lmax=0.png")

L2normerr = L2.L2_diff(V_cube, phi, g.scratch_array())
print("resolution", nx)
print("lmax =", n_moments)
print(L2normerr)
