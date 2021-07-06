"""
A simple model of the LSST focal plane
"""

import numpy as np
RADIAN_PER_DEGREE = np.pi / 180

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
        for ic in range(3):
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
