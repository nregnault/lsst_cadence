"""
A simple model of the LSST focal plane
"""

import numpy as np
from . import RADIAN_PER_DEGREE

class Chip(object):
    """
    LSST sensor. It knows 
       (1) where it is (in raft coordinates)
       (2) how to superpixellize itself. 

    .. example:
        >>> ccd = Chip(5, 1467., 7728.)
        >>> cells = ccd.pixellize(3,3)        
    """
    Nx, Ny = 4072, 4072
    
    def __init__(self, num, dx, dy):
        self.num = num
        self.dx, self.dy = dx, dy
        self.transform = np.matrix([[1., 0., self.dx], [0., 1., self.dy], [0., 0., 1.]])
        
    def pixellize(self, nx, ny):
        """Build and return an array of superpixels [in raft coordinates]

        Each superpixel is identified by the coordinates of its 4
        corners.  These corners are transformed from CCD coordinates
        to raft coordinates.

        Args:
          nx (int): number of cells along the x axis
          ny (int): number of cells along the y axis

        Returns:
          cells (recarray): an array of cells
        """
        ncells = nx*ny
        dx = np.floor(np.linspace(0, self.Nx, nx+1))
        dy = np.floor(np.linspace(0, self.Ny, ny+1))
        gx,gy = np.meshgrid(dx,dy)
        x0, y0 = gx[:-1,:-1].ravel(), gy[:-1,:-1].ravel()
        x1, y1 = gx[1:,1:].ravel(), gy[1:,1:].ravel()
        
        N = len(x0.ravel())        
        O = np.ones(N, dtype=np.float)
        p = np.zeros(N, dtype=[('iexp', 'i8'), ('icell', 'i4'), ('ichip', 'i4'), ('iraft', 'i4'), 
                               ('p0', '3f4'), ('p1', '3f4'), ('p2', '3f4'), ('p3', '3f4')])
        
        # initialize the cell structure 
        p['icell'] = np.arange(N, dtype=np.int)
        p['ichip'] = np.asarray(self.num * np.ones(N), dtype=np.int)
        p['iraft'] = -np.ones(N, dtype=np.int)
        p['p0'][:,0] = x0 ; p['p0'][:,1] = y0 ; p['p0'][:,2] = 1.
        p['p1'][:,0] = x1 ; p['p1'][:,1] = y0 ; p['p1'][:,2] = 1.
        p['p2'][:,0] = x1 ; p['p2'][:,1] = y1 ; p['p2'][:,2] = 1.
        p['p3'][:,0] = x0 ; p['p3'][:,1] = y1 ; p['p3'][:,2] = 1.
        
        # position it in the focal plane 
        M = self.transform
        for k in ['p0', 'p1', 'p2', 'p3']:
            p[k] = (M * p[k].T).T
        return p
        
    
class Raft(object):
    """
    A 3x3 ccd unit. 21 rafts on the LSST focal plane. 
    
    This class knows:
      (1) where it is w.r.t. the other rafts (in focal plane coordinates)
      (2) how to position its CCDs
      (3) how to superpixellize itself 
          (return its superpixels it focal plane coordinates).
    """
    Chip_gap = 100
    Nx, Ny = 3*Chip.Nx + 2*Chip_gap, 3*Chip.Ny + 2*Chip_gap
    
    CHIP_COORDINATES_IN_RAFT = np.asarray([
            (0, 0., 2.), (1, 1., 2.), (2, 2., 2.),
            (3, 0., 1.), (4, 1., 1.), (5, 2., 1.),
            (6, 0., 0.), (7, 1., 0.), (8, 2., 0.),
            ])
    
    def __init__(self, num, dx, dy):
        """
        """
        self.num = num
        self.dx, self.dy = dx, dy
        self.chips = []
        nnx, nny = Chip.Nx+self.Chip_gap, Chip.Ny+self.Chip_gap
        self.transform = np.matrix([[1., 0., self.dx], [0., 1., self.dy], [0., 0., 1.]])
        for ichip, dx, dy in self.CHIP_COORDINATES_IN_RAFT:
            self.chips.append(Chip(ichip, dx*nnx, dy*nny))
            
    def pixellize(self, nx, ny):
        """Superpixellize the raft. Superpixel Focal plane coordinates.
        
        Calls superpixellize on each of its CCDs, and transform the
        superpixel coordinates from Raft coordinates to CCD
        coordinates.
        """        
        p = np.hstack([c.pixellize(nx, ny) for c in self.chips])
        p['iraft'] = self.num
        # position the raft in the focal plane
        M = self.transform 
        for k in ['p0', 'p1', 'p2', 'p3']:
            p[k] = (M * p[k].T).T
        return p
    

class FocalPlane(object):
    """
    The full LSST focal plane (a set of 21 rafts)
    """    
    Raft_gap = 200
    N_rafts = 21
    N_ccd_per_raft = 9
    
    RAFT_COORDINATES_IN_FOCAL_PLANE = np.asarray([
            ( 0, -1.,  2.), ( 1,  0., 2.), ( 2, +1,  2.),
            ( 3, -2., 1.), ( 4, -1.,  1.), ( 5,  0,  1.), ( 6, +1., 1.), ( 7, +2., 1.),
            ( 8, -2., 0.), ( 9, -1.,  0.), (10,  0,  0.), (11, +1., 0.), (12, +2., 0.),
            (13, -2.,-1.), (14, -1., -1.), (15,  0, -1.), (16, +1.,-1.), (17, +2.,-1.),
            (18, -1., -2.), (19,  0.,-2.), (20, +1, -2.),])

    def __init__(self, nx=1, ny=1):
        self.pixel_size = 0.2
        self.nx, self.ny = 4072, 4072
        self.gap = 153
        self.rafts = []
        rad_per_degree = np.pi / 180.
        scale = 0.2 / 3600. * RADIAN_PER_DEGREE
        self.transform = np.matrix([[scale, 0., -6208.*scale], [0., scale, -6208. * scale], [0., 0., 1.]])
        nnx, nny = Raft.Nx + self.Raft_gap, Raft.Ny + self.Raft_gap
        for iraft, dx, dy in self.RAFT_COORDINATES_IN_FOCAL_PLANE:
            self.rafts.append(Raft(iraft, dx*nnx, dy*nny))
    
    def pixellize(self, nx, ny):
        p = np.hstack([r.pixellize(nx, ny) for r in self.rafts])
        # renumber the chips and the cells
        p['ichip'] += p['iraft'] * self.N_ccd_per_raft
        ncells_per_chip = nx*ny
        p['icell'] += p['ichip']*ncells_per_chip
        # recenter the mosaic, so that the center of raft 10
        # is at (0,0) in the tangent plane
        M = self.transform
        p['p0'] = (M * p['p0'].T).T
        p['p1'] = (M * p['p1'].T).T
        p['p2'] = (M * p['p2'].T).T
        p['p3'] = (M * p['p3'].T).T        
        return p

    def rotate(self, cells, theta=0.):
        t = theta * RADIAN_PER_DEGREE
        c,s = np.cos(t), np.sin(t)
        R = np.matrix([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
        cells['p0'] = (R * cells['p0'].T).T
        cells['p1'] = (R * cells['p1'].T).T
        cells['p2'] = (R * cells['p2'].T).T
        cells['p3'] = (R * cells['p3'].T).T
        return cells
        
    def tan2radec(self, xyz, tp):
        ra0, dec0  = (tp[0]  * RADIAN_PER_DEGREE, tp[1] * RADIAN_PER_DEGREE)
        cos0, sin0 = np.cos(dec0), np.sin(dec0)
        l,m,_ = xyz.T
        dec_t = cos0 - m * sin0
        ra_t = (ra0 + np.arctan2(l, dec_t) + np.pi) % (2. * np.pi) - np.pi                 
        dec_t = np.arctan(np.cos(ra_t-ra0) * (m*cos0 + sin0) / dec_t)        
        return ra_t, dec_t

    def to_sky(self, cells, tp):
        """
        takes as an input the cells defined in the TP and the (ra,dec)
        of the mosaic central point, and computes the (ra,dec) of the
        cell corners.
        
        .. note :: modifies the structure in place. 
        """
        for k in ['p0', 'p1', 'p2', 'p3']:
            ra, dec = self.tan2radec(cells[k], tp)
            cells[k][:,0] = ra
            cells[k][:,1] = dec
    
    def to_uv(self, cells, tp):
        """
        """
        for k in ['p0', 'p1', 'p2', 'p3']:
            ra, dec = self.tan2radec(cells[k], tp)
            cells[k][:,0] = np.cos(dec) * np.cos(ra)
            cells[k][:,1] = np.cos(dec) * np.sin(ra)
            cells[k][:,2] = np.sin(dec)


class FocalPlaneSimple(FocalPlane):
    """
    The LSST focal plane, with no sub-pixels.
    
    This class is used to study the LSST cadences (i.e. to derive from
    the cadence files, the effective cadence of each healpix pixel).
    
    """
    def __init__(self):
        """
        
        """
        p = self.cells = np.zeros(3,
                                  dtype=[('iexp', 'i8'), ('icell', 'i4'), ('ichip', 'i4'), ('iraft', 'i4'), 
                                         ('p0', '3f4'), ('p1', '3f4'), ('p2', '3f4'), ('p3', '3f4')])
        cells = {0: (-0.03048508, -0.01825226, -0.01825226, 0.01825226),
                 1: (-0.01825226,  0.01825226, -0.03048508, 0.03048508),
                 2: ( 0.01844619,  0.03048508, -0.01825226, 0.01825226)}
        self.rafts = []
        for ic in xrange(3):
            p[ic]['icell'] = ic
            xmin, xmax, ymin, ymax = cells[ic]
            p[ic]['p0'][:] = (xmin, ymin, 1.)
            p[ic]['p1'][:] = (xmax, ymin, 1.)
            p[ic]['p2'][:] = (xmax, ymax, 1.)
            p[ic]['p3'][:] = (xmin, ymax, 1.)        
    
    def pixellize(self):
        """
        Only one pixellization.
        """
        return self.cells.copy()
    

def draw_fp(cells, icell=None, newfig=True, alpha=0.5):
    import pylab as pl
    import matplotlib
    from matplotlib.patches import Rectangle, Polygon
    from matplotlib.collections import PatchCollection
    
    if newfig:
        fig = pl.figure(figsize=(8,8))
        sp = fig.add_subplot(111)        
    else:
        fig = pl.gcf()
        sp = pl.gca()
    l = []
    col = []
    for c in cells:
        if c['icell'] == icell:
            color = 1.
        else:
            color = c['icell']
        col.append(color)
        l.append(Polygon((c['p0'][:2], #* DEGREE_PER_RADIAN, 
                          c['p1'][:2], #* DEGREE_PER_RADIAN, 
                          c['p2'][:2], #* DEGREE_PER_RADIAN, 
                          c['p3'][:2], #* DEGREE_PER_RADIAN
                          )))
    pc = PatchCollection(l, alpha=0.5, lw=2)
    sp.add_collection(pc)
    pc.set_array(np.array(col))
    pl.xlim((cells['p0'][:,0].min()-0.01, cells['p2'][:,0].max()+0.01))
    pl.ylim((cells['p0'][:,1].min()-0.01, cells['p2'][:,1].max()+0.01))
    pl.show()
