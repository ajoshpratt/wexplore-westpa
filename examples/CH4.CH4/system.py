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
    ncenters = centers.shape[0]
    d = np.empty((ncenters,), dtype=np.float32)

    for k in xrange(ncenters):
        d[k] = np.sqrt((p[0] - centers[k,0])**2 + (p[1] - centers[k,1])**2)

    return d

class System(WESTSystem):
    def initialize(self):
        self.pcoord_ndim = 2
        self.pcoord_len = 11
        self.pcoord_dtype = numpy.float32
        self.dist_binbounds = [10.00,float('inf')]
        self.rmsd_bounds = [0.0, float('inf')]
        self.wexploreMapper = wexplore.WExploreBinMapper(n_regions=[4,4,4], d_cut=[2, .75, 0.25], dfunc=dfunc)
        self.wexploreMapper.centers = [[2.5, 7.0]]
        self.wexploreMapper.add_bin(None, 0)
        self.wexploreMapper.mask = np.array([[True, False]])
        self.wexploreMapper.bin_target_counts = 64
        self.initialize_mappers(self.wexploreMapper)
        self.max_replicas = 64
        self.bin_mapper.centers = [[2.5, 7.0]]
        self.tstates = ['unbound,10.0,1.0']
        self.n_regions = [6,4,4]

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

