import numpy as np


class Grid():
    def __init__(self, nx, ny, nz, xlim, ylim, zlim):
        # xlim=(0.0, 1.0), ylim=(0.0, 1.0), zlim=(-1.0, 1.0)
        """an axisymmetric grid"""

        self.nx = nx
        self.ny = ny
        self.nz = nz

        self.xlim = xlim
        self.ylim = ylim
        self.zlim = zlim

        self.dx = (self.xlim[1] - self.xlim[0])/self.nx
        self.dy = (self.ylim[1] - self.ylim[0])/self.ny
        self.dz = (self.zlim[1] - self.zlim[0])/self.nz

        """
        self.nx = nx-1
        self.ny = ny-1
        self.nz = nz-1

        self.xlim = xlim
        self.ylim = ylim
        self.zlim = zlim

        self.dx = (self.xlim[1] - self.xlim[0])/nx
        self.dy = (self.ylim[1] - self.ylim[0])/ny
        self.dz = (self.zlim[1] - self.zlim[0])/nz
        """

        self.x = (np.arange(self.nx) + 0.5)*self.dx + self.xlim[0]
        self.y = (np.arange(self.ny) + 0.5)*self.dy + self.ylim[0]
        self.z = (np.arange(self.nz) + 0.5)*self.dz + self.zlim[0]

        self.x3d, self.y3d, self.z3d = np.meshgrid(self.x, self.y, self.z,
                                                   indexing='ij')

        self.dr = np.sqrt(self.dx**2+self.dy**2+self.dz**2)

        # cartisian volume element
        self.vol = self.dx*self.dy*self.dz

    def scratch_array(self):
        return np.zeros((self.nx, self.ny, self.nz), dtype=np.float64)

    def scratch_x_plane_array(self):
        return np.zeros((self.ny, self.nz), dtype=np.float64)

    def scratch_y_plane_array(self):
        return np.zeros((self.nx, self.nz), dtype=np.float64)

    def scratch_z_plane_array(self):
        return np.zeros((self.nx, self.ny), dtype=np.float64)
