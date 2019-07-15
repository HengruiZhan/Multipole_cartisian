import numpy as np
from scipy.special import sph_harm


class Multipole():
    def __init__(self, grid, rho, l_moments, dr, center=(0.0, 0.0)):

        self.g = grid
        self.l_moments = l_moments
        self.dr_mp = dr
        self.center = center

        # compute the bins
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
        # we'll index the list by multipole moment l
        self.m_r = []
        self.m_i = []

        for l in range(l_moments):
            self.m_r.append([])
            self.m_i.append([])
        for l in range(l_moments):
            for m in range(2*l+1):
                self.m_r[l].append(np.zeros(self.n_bins, dtype=np.complex128))
                self.m_i[l].append(np.zeros(self.n_bins, dtype=np.complex128))
        for l in range(l_moments):
            for m in range(-l, l+1):
                self.compute_expansion(rho, l, m)

    def compute_expansion(self, rho, l, m):
        # rho is density that lives on a grid self.g

        # loop over cells
        for i in range(self.g.nx):
            for j in range(self.g.ny):
                for k in range(self.g.nz):

                    # for each cell, i,j, compute r, theta (polar angle
                    # from z) and phi(azumuthal angle)
                    # and determine which shell we are in
                    radius = np.sqrt((self.g.x[i] - self.center[0])**2 +
                                     (self.g.y[j] - self.center[1])**2 +
                                     (self.g.z[k] - self.center[2])**2)
                    r = np.sqrt((self.g.x[i] - self.center[0])**2 +
                                (self.g.y[j] - self.center[1])**2)
                    # tan(theta) = r/z
                    theta = np.arctan2(r, self.g.z[k])

                    # tan(phi) = y/x
                    phi = np.arctan2(self.g.y[j], self.g.x[i])

                    # loop over the multipole moments, l (m = 0 here)
                    m_zone = rho[i, j, k] * self.g.vol

                    # compute Y_l^m (note: we use theta as the polar
                    # angle, phi as the azumuthal angle, scipy is opposite)
                    Y_lm = sph_harm(m, l, phi, theta)

                    R_lm = np.sqrt(4*np.pi/(2*l + 1)) * radius**l * Y_lm
                    I_lm = np.sqrt(4*np.pi/(2*l + 1)) * Y_lm / radius**(l+1)

                    # add to the all of the appropriate inner or outer
                    # moment functions
                    imask = radius <= self.r_bin
                    omask = radius > self.r_bin

                    # 0 <= l, -l <= m <= l
                    # -l+l=0 corresponds to m= -l, l+l=2l corresponds to +l
                    self.m_r[l][m+l][imask] += R_lm * m_zone
                    self.m_i[l][m+l][omask] += I_lm * m_zone

    def sample_mtilde(self, r, l, m):
        # this returns the result of Eq. 19

        # we need to find which be we are in
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

    def phi(self, x, y, z):
        # return Phi(r), using Eq. 20

        radius = np.sqrt((x - self.center[0])**2 +
                         (y - self.center[1])**2 +
                         (z - self.center[2])**2)

        r = np.sqrt((x - self.center[0])**2 +
                    (y - self.center[1])**2)

        # tan(theta) = r/z
        theta = np.arctan2(r, z)

        # tan(phi) = y/x
        phi = np.arctan2(y, x)

        phi_zone = 0.0

        for l in range(self.l_moments):
            for m in range(-l, l+1):
                mtilde_r, mtilde_i = self.sample_mtilde(radius, l, m)

                Y_lm = sph_harm(m, l, phi, theta)
                R_lm = np.sqrt(4*np.pi/(2*l + 1)) * radius**l * Y_lm
                I_lm = np.sqrt(4*np.pi/(2*l + 1)) * Y_lm / radius**(l+1)

                phi_zone += mtilde_r * np.conj(I_lm) + np.conj(mtilde_i) * R_lm

        return -np.real(phi_zone)


"""
class Multipole():
    def __init__(self, grid, l_moments, dr, center=(0.0, 0.0)):

        self.g = grid
        self.l_moments = l_moments
        self.dr_mp = dr
        self.center = center

        # compute the bins
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
        # we'll index the list by multipole moment l
        self.m_r = []
        self.m_i = []

        for l in range(self.l_moments):
            self.m_r.append(np.zeros((self.n_bins), dtype=np.complex128))
            self.m_i.append(np.zeros((self.n_bins), dtype=np.complex128))

    # def compute_expansion(self, rho, l, m):
    def compute_expansion(self, rho, l):
        # rho is density that lives on a grid self.g

        # loop over cells
        for i in range(self.g.nx):
            for j in range(self.g.ny):
                for k in range(self.g.nz):

                    # for each cell, i,j, compute r, theta (polar angle
                    # from z) and phi(azumuthal angle)
                    # and determine which shell we are in
                    radius = np.sqrt((self.g.x[i] - self.center[0])**2 +
                                     (self.g.y[j] - self.center[1])**2 +
                                     (self.g.z[k] - self.center[2])**2)

                    r = np.sqrt((self.g.x[i] - self.center[0])**2 +
                                (self.g.y[j] - self.center[1])**2)
                    # tan(theta) = r/z
                    theta = np.arctan2(r, self.g.z[k])

                    # tan(phi) = y/x
                    # phi = np.arctan2(self.g.y[j], self.g.x[i])

                    # loop over the multipole moments, l (m = 0 here)
                    m_zone = rho[i, j, k] * self.g.vol

                    # compute Y_l^m (note: we use theta as the polar
                    # angle, phi as the azumuthal angle, scipy is opposite)
                    # Y_lm = sph_harm(m, l, phi, theta)
                    Y_lm = sph_harm(0, l, 0.0, theta)

                    R_lm = np.sqrt(4*np.pi/(2*l + 1)) * radius**l * Y_lm
                    I_lm = np.sqrt(4*np.pi/(2*l + 1)) * Y_lm / radius**(l+1)

                    # add to the all of the appropriate inner or outer
                    # moment functions
                    imask = radius <= self.r_bin
                    omask = radius > self.r_bin

                    self.m_r[l][imask] += R_lm * m_zone
                    self.m_i[l][omask] += I_lm * m_zone

    # def sample_mtilde(self, r, l, m):
    def sample_mtilde(self, l, r):
        # this returns the result of Eq. 19

        # we need to find which be we are in
        mu_m = np.argwhere(self.r_bin <= r)[-1][0]
        mu_p = np.argwhere(self.r_bin > r)[0][0]

        assert mu_p == mu_m + 1

        mtilde_r = (r - self.r_bin[mu_m])/(self.r_bin[mu_p] - self.r_bin[mu_m]
                                           ) * self.m_r[l][mu_p] + \
            (r - self.r_bin[mu_p])/(self.r_bin[mu_m] -
                                    self.r_bin[mu_p]) * self.m_r[l][mu_m]

        mtilde_i = (r - self.r_bin[mu_m])/(self.r_bin[mu_p] - self.r_bin[mu_m]
                                           ) * self.m_i[l][mu_p] + \
            (r - self.r_bin[mu_p])/(self.r_bin[mu_m] -
                                    self.r_bin[mu_p]) * self.m_i[l][mu_m]

        return mtilde_r, mtilde_i

    def phi(self, x, y, z):
        # return Phi(r), using Eq. 20

        radius = np.sqrt((x - self.center[0])**2 +
                         (y - self.center[1])**2 +
                         (z - self.center[2])**2)

        r = np.sqrt((x - self.center[0])**2 +
                    (y - self.center[1])**2)

        # tan(theta) = r/z
        theta = np.arctan2(r, z)

        # tan(phi) = y/x
        # phi = np.arctan2(y, x)

        phi_zone = 0.0

        for l in range(self.l_moments):
            mtilde_r, mtilde_i = self.sample_mtilde(l, radius)

            Y_lm = sph_harm(0, l, 0.0, theta)
            R_lm = np.sqrt(4*np.pi/(2*l + 1)) * radius**l * Y_lm
            I_lm = np.sqrt(4*np.pi/(2*l + 1)) * Y_lm / radius**(l+1)

            phi_zone += mtilde_r * np.conj(I_lm) + np.conj(mtilde_i) * R_lm

        return -np.real(phi_zone)
        """
