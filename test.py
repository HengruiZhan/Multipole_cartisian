import analytical_potential as ap
import grid
import matplotlib.pyplot as plt
import multipole
import numpy as np
import L2_difference as L2
# from gekko import GEKKO
"""
nx = 128
ny = 128
nz = 128
"""
nx = 100
ny = 100
nz = 100
xlim = (-1.0, 1.0)
ylim = (-1.0, 1.0)
zlim = (-1.0, 1.0)
g = grid.Grid(nx, ny, nz, xlim, ylim, zlim)
dens = g.scratch_array()
"""
# density of the sphere
sph_center = (0.0, 0.0, 0.3)
radius = np.sqrt((g.x3d - sph_center[0])**2 + (g.y3d - sph_center[1])**2 +
                 (g.z3d - sph_center[2])**2)
dens[radius <= 0.3] = 1.0
"""

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
# plot the density
plt.imshow(np.transpose(dens), origin="lower",
           interpolation="nearest",
           extent=[g.rlim[0], g.rlim[1],
                   g.zlim[0], g.zlim[1]])

ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("dens.png")
"""

n_moments = 4
# center = (a/2, b/2, c/2)
center = (0.0, 0.0, 0.0)
moment = multipole.Multipole(g, dens, n_moments, 2*g.dr, center=center)

"""
# with azumuthal
for l in range(l_moments):
    for m in range(-l, l+1):
        moment.compute_expansion(dens, l, m)
"""

phi = g.scratch_y_plane_array()

for i in range(g.nx):
    for j in range(g.nz):
        # phi[i, j] = moment.phi(g.x[i] - g.dx/2, g.y[0], g.z[j])
        phi[i, j] = moment.phi(g.x[i], g.y[0], g.z[j], g.dx/2)

plt.imshow(np.log10(np.abs(np.transpose(phi))), origin="lower",
           interpolation="nearest",
           extent=[g.xlim[0], g.xlim[1],
                   g.zlim[0], g.zlim[1]])

plt.colorbar()
ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("phi.png")

L2normerr = L2.L2_diff(V_cube, phi, g.scratch_array())
print("x_off with dx/2")
print("lmax =", n_moments-1)
print(L2normerr)

"""
xp = np.array([0, 1, 2, 3, 4, 5, 6])
yp = np.array([0.004129259053218379, 0.003912801608538426, 0.00391070661969043,
               0.003913627885664768, 0.003913372974168947,
               0.003913363646001866, 0.003913366211304788])
x = np.arange(0, 22)
y = np.interp(x, xp, yp)
plt.plot(xp, yp, 'o')
plt.plot(x, y, 'r--')
plt.show()
"""
