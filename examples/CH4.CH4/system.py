from __future__ import division, print_function; __metaclass__ = type
import os, sys, math, itertools
import numpy
import numpy as np
import west
from west import WESTSystem
from westpa.binning import RectilinearBinMapper, RecursiveBinMapper, FuncBinMapper
import westpa

sys.path.append(os.path.dirname(__file__) + "wexplore")
from wexplore import wexplore, wex_utils
from scipy.spatial.distance import cdist

import logging
log = logging.getLogger(__name__)
log.debug('loading module %r' % __name__)

coord_dtype = numpy.float32

def dfunc_normal(coordvec, centers):
    if coordvec.ndim < 2:
        new_coordvec = np.empty((1,coordvec.shape[0]), dtype=coord_dtype)
        new_coordvec[0,:] = coordvec[:]
        coordvec = new_coordvec 
    distmat = np.require(cdist(coordvec, centers), dtype=coord_dtype)

    return distmat[0,:]

def fn2(coords, mask, output):
    test = coords[:,0] < 10.0
    # These are the ones for which test is true.
    output[mask & test] = 0
    output[mask & ~test] = 1

def dfunc(p, centers):
    # Should be a function which always returns the distance in the second dimension as 0.
    # This is useful for 'cancelling' out a dimension, if necessary.
    ncenters = centers.shape[0]
    d = np.empty((ncenters,), dtype=np.float32)

    for k in xrange(ncenters):
        d[k] = np.sqrt((p[0] - centers[k,0])**2 + (p[1] - centers[k,1])**2)
        #d[k] = np.sqrt((p[0] - centers[k,0])**2)
    #d[1] = 0

    return d

def dfunc_filter(p, centers):
    # Should be a function which always returns the distance in the second dimension as 0.
    # This is useful for 'cancelling' out a dimension, if necessary.
    ncenters = centers.shape[0]
    d = np.empty((ncenters,), dtype=np.float32)

    if p[1] <= 5.0:
        for k in xrange(ncenters):
            d[k] = np.sqrt((p[0] - centers[k,0])**2 + (p[1] - centers[k,1])**2)
            #d[k] = np.sqrt((p[0] - centers[k,0])**2)
        #d[1] = 0
    else:
        for k in xrange(ncenters):
            d[k] = np.sqrt((p[1] - centers[k,1])**2)

    return d

class System(WESTSystem):
    def initialize(self):
        self.pcoord_ndim = 2
        self.pcoord_len = 11
        self.pcoord_dtype = numpy.float32
        #self.dist_binbounds = [0.0,2.80,2.88,3.00,3.10,3.29,3.79,3.94,4.12,4.39,5.43,5.90,6.90,7.90,8.90,9.90,10.90,11.90,12.90,13.90,14.98,15.90,float('inf')]
        #self.dist_binbounds = [7.00,7.75,8.50,9.25,10.00,float('inf')]
        self.dist_binbounds = [10.00,float('inf')]
        self.rmsd_bounds = [0.0, float('inf')]
        #self.color_binbounds = [0,1,2,float('inf')]
        self.wexploreMapper = wexplore.WExploreBinMapper(n_regions=[4,4,4], d_cut=[2, .75, 0.25], dfunc=dfunc)
        self.wexploreMapper.centers = [[2.5, 7.0]]
        self.wexploreMapper.add_bin(None, 0)
        self.wexploreMapper.mask = np.array([[True, False]])
        self.wexploreMapper.bin_target_counts = 64
        self.initialize_mappers(self.wexploreMapper)
        self.max_replicas = 64
        self.bin_mapper.centers = [[2.5, 7.0]]
        self.tstates = ['unbound,10.0,1.0']

    def initialize_mappers(self,wexploreMapper):
        self.outerMapper = FuncBinMapper(fn2,2)
        self.recurseMapper = RecursiveBinMapper(self.outerMapper)
        self.RectiMapper = RectilinearBinMapper([self.dist_binbounds, self.rmsd_bounds])
        self.RectiMapper.bin_target_counts = 3
        self.recurseMapper.add_mapper(wexploreMapper, [4.99, 0])
        self.recurseMapper.add_mapper(self.RectiMapper, [10, 0])
        self.bin_mapper = self.recurseMapper
        self.bin_target_counts = numpy.empty((self.bin_mapper.nbins,), numpy.int)
        self.bin_target_counts[...] = 64
        return self.bin_mapper

    
def coord_loader(fieldname, coord_file, segment, single_point=False):
    coord_raw = numpy.loadtxt(coord_file, dtype=numpy.float32) 
    color_bins = [(0.0,4.00),(10.00,float('inf'))]
    unknown_state = 2
    new_walker = 0
    system = westpa.rc.get_system_driver()
    #print(system.pcoord_ndim)

    try:
        npts = len(coord_raw)
    except(TypeError):
        npts = 1
        new_walker = 1

    coords = numpy.empty((npts), numpy.float32)
    colors = numpy.empty((npts), numpy.float32)
    if new_walker == 1:
        colors[:] = unknown_state
        for istate,state_tuple in enumerate(color_bins):
            if coord_raw >= state_tuple[0] and coord_raw < state_tuple[1]:
                colors[:] = istate
    else:
        colors[:] = segment.pcoord[0][1]
    if new_walker == 1:
        coords[0] = coord_raw
    else:
        coords[:] = coord_raw
    for istate,state_tuple in enumerate(color_bins):
        if coords[-1] >= state_tuple[0] and coords[-1] < state_tuple[1]:
            colors[-1] = istate
    
    if new_walker == 1:
        segment.pcoord = numpy.hstack((coords[:],colors[:]))
    else:
        segment.pcoord = numpy.swapaxes(numpy.vstack((coords[:],colors[:])), 0, 1)

    
    
    
    
    
