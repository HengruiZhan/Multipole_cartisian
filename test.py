import analytical_potential as ap
import grid
import matplotlib.pyplot as plt
import multipole
import numpy as np

nx = 128
ny = 128
nz = 256
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
V_cube = g.scratch_y_plane_array()
for i in range(g.nx):
    for j in range(g.nz):
        cube_V = ap.Analytical_Cube_Potential(sites, g.x[i], g.y[0], g.z[j])
        V_cube[i, j] = cube_V.V
plt.imshow(np.log10(np.abs(np.transpose(V_cube))), origin="lower",
           interpolation="nearest",
           extent=[g.xlim[0], g.xlim[1],
                   g.zlim[0], g.zlim[1]])

plt.colorbar()
ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("phi.png")
"""
for i in range(g.nx):
    for j in range(g.ny):
        for k in range(g.nz):
            if g.x[i] >= 0 and \
                    g.x[i] <= a and \
                    g.y[j] >= 0 and \
                    g.y[j] <= b and \
                    g.z[k] >= 0 and \
                    g.z[k] <= c:
                dens[i, j, k] = 1
"""
"""
plt.imshow(np.transpose(dens), origin="lower",
           interpolation="nearest",
           extent=[g.rlim[0], g.rlim[1],
                   g.zlim[0], g.zlim[1]])

ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("dens.png")
"""
# without azumuthal
# l_moments = 3
# center = (a/2, b/2, c/2)
# center = (0.0, 0.0, 0.0)
# moment = multipole.Multipole(g, 3, 2*g.dr, center=center)
# moment = multipole.Multipole(g, dens, 3, 2*g.dr, center=center)
# with azumuthal
"""
for l in range(l_moments):
    for m in range(-l, l+1):
        moment.compute_expansion(dens, l, m)
"""
"""
for l in range(l_moments):
    moment.compute_expansion(dens, l)
    """
# m.compute_expansion(dens, 1, -1)

# phi = g.scratch_y_plane_array()
"""
for i in range(g.nx):
    for j in range(g.nz):
        phi[i, j] = moment.phi(g.x[i], g.y[0], g.z[j])

plt.imshow(np.log10(np.abs(np.transpose(phi))), origin="lower",
           interpolation="nearest",
           extent=[g.xlim[0], g.xlim[1],
                   g.zlim[0], g.zlim[1]])

plt.colorbar()
ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("phi.png")
"""
