# newer version of multipole for 2d axisymetric case
from numba import jit
import numpy as np
import scipy.constants as sc
from scipy.special import sph_harm


class Multipole():
    def __init__(self, grid, density, l_moments, dr, center=(0.0, 0.0)):
        # multipole approximation of the potential

        self.g = grid
        self.l_moments = l_moments+1
        self.dr_mp = dr
        self.center = center

        # compute the bins, or the radius of the concentric sphere, r_mu
        x_max = max(abs(self.g.xlim[0] - center[0]), abs(self.g.xlim[1] -
                                                         center[0]))
        y_max = max(abs(self.g.ylim[0] - center[0]), abs(self.g.ylim[1] -
                                                         center[1]))
        z_max = max(abs(self.g.zlim[0] - center[1]), abs(self.g.zlim[1] -
                                                         center[2]))

        dmax = np.sqrt(x_max**2 + y_max**2 + z_max**2)

        self.n_bins = int(dmax/dr)

        # bin boundaries
        self.r_bin = np.linspace(0.0, dmax, self.n_bins)

        # storage for the inner and outer multipole moment functions
        # we'll index the list by multipole moment l and m
        self.m_r = []
        self.m_i = []

        self.m_r = []
        self.m_i = []
        for l in range(self.l_moments):
            self.m_r.append(np.zeros((self.n_bins), dtype=np.complex128))
            self.m_i.append(np.zeros((self.n_bins), dtype=np.complex128))
        for l in range(l_moments):
            self.compute_expansion(density, l, 0)
        """
        # modify and only consider l, to simplify the calculation
        for l in range(l_moments):
            self.m_r.append([])
            self.m_i.append([])
        # for each l, there are 2l+1 values for m, -l <= m <= +l
        for l in range(l_moments):
            for m in range(2*l+1):
                self.m_r[l].append(np.zeros(self.n_bins, dtype=np.complex128))
                self.m_i[l].append(np.zeros(self.n_bins, dtype=np.complex128))
        # calculate the m_r and m_i for every l and m
        for l in range(l_moments):
            for m in range(-l, l+1):
                self.compute_expansion(density, l, m)
                """

    @jit
    def compute_expansion(self, density, l, m):
        # this calculate the m_r and m_i indexed l, m,
        # defined in eq. 17 and eq. 18
        # density is density that lives on a grid self.g

        # loop over cells
        for i in range(self.g.nx):
            for j in range(self.g.ny):
                for k in range(self.g.nz):

                    # for each cell, i,j,k compute radius, rho, theta (polar
                    # angle) and phi(azumuthal angle)
                    # and determine which shell we are in
                    radius = np.sqrt((self.g.x[i] - self.center[0])**2 +
                                     (self.g.y[j] - self.center[1])**2 +
                                     (self.g.z[k] - self.center[2])**2)
                    rho = np.sqrt((self.g.x[i] - self.center[0])**2 +
                                  (self.g.y[j] - self.center[1])**2)
                    # tan(theta) = r/z
                    theta = np.arctan2(rho, self.g.z[k])

                    # tan(phi) = y/x
                    phi = np.arctan2(self.g.y[j], self.g.x[i])

                    # loop over the multipole moments, l (m = 0 here)
                    m_zone = density[i, j, k] * self.g.vol

                    # compute Y_l^m (note: we use theta as the polar
                    # angle, phi as the azumuthal angle, scipy is opposite)
                    Y_lm = sph_harm(m, l, phi, theta)

                    R_lm = np.sqrt(4*np.pi/(2*l + 1)) * radius**l * Y_lm
                    I_lm = np.sqrt(4*np.pi/(2*l + 1)) * Y_lm / radius**(l+1)

                    # add to the all of the appropriate inner or outer
                    # moment functions
                    imask = radius <= self.r_bin
                    omask = radius > self.r_bin

                    """
                    # 0 <= l, -l <= m <= l
                    # -l+l=0 corresponds to m=-l, l+l=2l corresponds to m=+l
                    self.m_r[l][m+l][imask] += R_lm * m_zone
                    self.m_i[l][m+l][omask] += I_lm * m_zone
                    """
                    self.m_r[l][imask] += R_lm * m_zone
                    self.m_i[l][omask] += I_lm * m_zone

    def sample_mtilde(self, r, l, m):
        # this returns the result of Eq. 19
        # r is the radius of the point of the field from the expansion center

        # we need to find out the index of r_mu_plus and r_mu_minus in Eq. 19
        mu_m = np.argwhere(self.r_bin <= r)[-1][0]
        mu_p = np.argwhere(self.r_bin > r)[0][0]

        assert mu_p == mu_m + 1

        mtilde_r = (r - self.r_bin[mu_m])/(self.r_bin[mu_p] - self.r_bin[mu_m]
                                           ) * self.m_r[l][m+l][mu_p] + \
            (r - self.r_bin[mu_p])/(self.r_bin[mu_m] -
                                    self.r_bin[mu_p]) * self.m_r[l][m+l][mu_m]

        mtilde_i = (r - self.r_bin[mu_m])/(self.r_bin[mu_p] - self.r_bin[mu_m]
                                           ) * self.m_i[l][m+l][mu_p] + \
            (r - self.r_bin[mu_p])/(self.r_bin[mu_m] -
                                    self.r_bin[mu_p]) * self.m_i[l][m+l][mu_m]

        return mtilde_r, mtilde_i

    @jit
    def compute_phi(self, x, y, z, dx, dy, dz):
        # return Phi(r), potential of one point of the field, using Eq. 20
        # x, y, z are x-, y-, z- coordinates of the point of the field

        x -= dx
        y -= dy
        z -= dz

        radius = np.sqrt((x - self.center[0])**2 +
                         (y - self.center[1])**2 +
                         (z - self.center[2])**2)

        rho = np.sqrt((x - self.center[0])**2 +
                      (y - self.center[1])**2)

        # tan(theta) = r/z
        theta = np.arctan2(rho, z)

        # tan(phi) = y/x
        phi = np.arctan2(y, x)

        phi_zone = 0.0

        """
        for l in range(self.l_moments):
            for m in range(-l, l+1):
                mtilde_r, mtilde_i = self.sample_mtilde(radius, l, m)

                Y_lm = sph_harm(m, l, phi, theta)
                R_lm = np.sqrt(4*np.pi/(2*l + 1)) * radius**l * Y_lm
                I_lm = np.sqrt(4*np.pi/(2*l + 1)) * Y_lm / radius**(l+1)

                phi_zone += mtilde_r * np.conj(I_lm) + np.conj(mtilde_i) * R_lm
                """
        # only consider l when axysymetric
        for l in range(self.l_moments):
            mtilde_r, mtilde_i = self.sample_mtilde(radius, l, 0)

            Y_lm = sph_harm(0, l, phi, theta)
            R_lm = np.sqrt(4*np.pi/(2*l + 1)) * radius**l * Y_lm
            I_lm = np.sqrt(4*np.pi/(2*l + 1)) * Y_lm / radius**(l+1)

            # phi_zone += mtilde_r * np.conj(I_lm) + np.conj(mtilde_i) * R_lm
            phi_zone += sc.G * (mtilde_r * np.conj(I_lm) +
                                np.conj(mtilde_i) * R_lm)

        return -np.real(phi_zone)

    def phi_point(self, x, y, z):
        """calculate the potential of a specific point"""
        dx = self.g.dx/2
        dy = self.g.dy/2
        dz = self.g.dz/2

        # phi_0 = self.compute_phi(r, z, 0, 0)
        phi_m_x = self.compute_phi(x, y, z, -dx, 0, 0)
        phi_m_y = self.compute_phi(x, y, z, 0, -dy, 0)
        phi_m_z = self.compute_phi(x, y, z, 0, 0, -dz)
        phi_p_x = self.compute_phi(x, y, z, dx, 0, 0)
        phi_p_y = self.compute_phi(x, y, z, 0, dy, 0)
        phi_p_z = self.compute_phi(x, y, z, 0, 0, dz)

        # phi = (phi_m_r+phi_m_z+phi_p_r+phi_p_z)/6

        phi = (phi_m_x+phi_m_y+phi_m_z+phi_p_x+phi_p_y+phi_p_z)/6

        return phi
