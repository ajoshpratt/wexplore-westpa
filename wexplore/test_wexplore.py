import numpy as np

from scipy.spatial.distance import cdist

import nose
import nose.tools

import westpa
import westtools
from wexplore import WExploreBinMapper
from westpa.binning import RectilinearBinMapper, RecursiveBinMapper, FuncBinMapper

coord_dtype = np.float32


def dfunc(coordvec, centers):
    if coordvec.ndim < 2:
        new_coordvec = np.empty((1,coordvec.shape[0]), dtype=coord_dtype)
        new_coordvec[0,:] = coordvec[:]
        coordvec = new_coordvec 
    distmat = np.require(cdist(coordvec, centers), dtype=coord_dtype)

    return distmat[0,:]

def fn2(coords, mask, output):
    # Arbitrary function for functional mapper
    test = coords[:,0] < 10.0
    output[mask & test] = 0
    output[mask & ~test] = 1

class TestWExploreBinMapper:

    def test_create_mapper(self):
        bin_mapper = WExploreBinMapper([2, 3, 3], [4.0, 2.0, 1.0], dfunc)

        assert bin_mapper.n_levels == 3
        assert bin_mapper.nbins == 0
        assert bin_mapper.max_nbins == 18

    def test_nested_mapper(self):
        # Tests for nested mapping.  Setting up the system (some arbitrary values, here)
        wexploreMapper = WExploreBinMapper(n_regions=[4,4,4], d_cut=[2, .75, 0.25], dfunc=dfunc)
        wexploreMapper.centers = [[2.5, 7.0]]
        wexploreMapper.add_bin(None, 0)

        outerMapper = FuncBinMapper(fn2,2)
        recurseMapper = RecursiveBinMapper(outerMapper)
        RectiMapper = RectilinearBinMapper([[10.00, float('inf')], [0.0, float('inf')]])
        recurseMapper.add_mapper(wexploreMapper, [9.9, 0])
        recurseMapper.add_mapper(RectiMapper, [10, 0])
        bin_mapper = recurseMapper
        
        # Now, begin testing setup...
        assert bin_mapper.nbins == 2

        wexploreMappers = []
        for key,i in bin_mapper.mapper_list.iteritems():
            if isinstance(i['base_mapper'], WExploreBinMapper) == True:
                wexploreMappers.append(i)
        self.wexploreMappers = wexploreMappers
        wexploreMappers[0]['base_mapper'].add_bin(0,1)

        bin_mapper.refresh_mappers()

        assert bin_mapper.nbins == 3

    def test_add_bin(self):
        bin_mapper = WExploreBinMapper([2, 2, 2], [4.0, 2.0, 1.0], dfunc)

        bin_mapper.add_bin(None, [0.0, 0.0])

        assert bin_mapper.nbins == 1
        assert bin_mapper.next_bin_index == 3
        for k in xrange(3):
            assert len(bin_mapper.level_indices[k]) == 1

        bin_mapper.add_bin(None, [-3.0, 0.0])

        assert bin_mapper.nbins == 2
        for k in xrange(3):
            assert len(bin_mapper.level_indices[k]) == 2

        node = bin_mapper.level_indices[1][0]
        bin_mapper.add_bin(node, [-0.5, 0.0])

        assert bin_mapper.nbins == 3
        assert len(bin_mapper.level_indices[0]) == 2
        assert len(bin_mapper.level_indices[1]) == 2
        assert len(bin_mapper.level_indices[2]) == 3

        # Should not add bins when adding a bin with a parent
        # in the lowest level
        node = bin_mapper.level_indices[2][0]
        bin_mapper.add_bin(node, [-0.5, 0.0])

        assert bin_mapper.nbins == 3
        assert len(bin_mapper.level_indices[0]) == 2
        assert len(bin_mapper.level_indices[1]) == 2
        assert len(bin_mapper.level_indices[2]) == 3

    def test_fetch_centers(self):
        bin_mapper = WExploreBinMapper([3, 3, 3], [4.0, 2.0, 1.0], dfunc)

        coords = np.arange(20).reshape((10,2))

        bin_mapper.centers = coords[:3]

        for k in xrange(3):
            bin_mapper.add_bin(None, k)

        centers = bin_mapper.fetch_centers(bin_mapper.level_indices[0])
        assert np.allclose(centers, coords[:3,:])

        centers = bin_mapper.fetch_centers([0, 1, 2])
        assert np.allclose(centers, 
                np.vstack((coords[0], coords[0], coords[0])))

        centers = bin_mapper.fetch_centers([0])
        assert np.allclose(centers, coords[0])

    def test_balance_replicas(self):
        bin_mapper = WExploreBinMapper([2, 2, 2], [4.0, 2.0, 1.0], dfunc)

        for k in xrange(2):
            bin_mapper.add_bin(None, None)

        for node in bin_mapper.level_indices[0]:
            bin_mapper.add_bin(node, None)

        for node in bin_mapper.level_indices[1]:
            bin_mapper.add_bin(node, None)

        assignments = np.array([1, 0, 0, 0, 0, 0, 0, 2, 2, 3, 4, 5, 6, 6, 6, 7])
        target_count = bin_mapper.balance_replicas(16, assignments)

        assert (target_count == 2*np.ones(8, dtype=np.int)).all()

    def test_assign_1D(self):
        bin_mapper = WExploreBinMapper([2, 2, 2], [4.0, 2.0, 1.0], dfunc)
        bin_mapper.centers = [[-1.0], [1.0], [-2.0], [2.0], [10.0]]

        bin_mapper.add_bin(None, 0)
        bin_mapper.add_bin(None, 1)

        for node in bin_mapper.level_indices[0]:
            cix = bin_mapper.bin_graph.node[node]['center_ix']
            if bin_mapper.centers[cix][0] > 0:
                bin_mapper.add_bin(node, 3)
            else:
                bin_mapper.add_bin(node, 2)

        for nix in [1, 4, 6, 8]:
            bin_mapper.add_bin(nix, 4)

        pcoords = np.array([[-1.0], [1.0], [2.0], 
                            [-2.0], [10.0], [-10.0], 
                            [3.0], [-3.0]], np.float32)
        assign = bin_mapper.assign(pcoords)
        assert list(assign) == [0, 1, 3, 2, 7, 2, 3, 2]

        assign = bin_mapper.assign(pcoords)
        assert list(assign) == [0, 1, 3, 2, 7, 2, 3, 2]

    def test_add_bin_distance(self):
        bin_mapper = WExploreBinMapper([2, 2, 2], [4.0, 2.0, 1.0], dfunc)
        bin_mapper.centers = [[0.0]]
        bin_mapper.add_bin(None, 0)

        pcoords = np.array([[0.0], [4.1], [5.1]], np.float32)
        assign = bin_mapper.assign(pcoords, add_bins=True)
        assert list(assign) == [0, 0, 0]
        assert bin_mapper.nbins == 2

        assign = bin_mapper.assign(pcoords)
        assert list(assign) == [0, 1, 1]

        #---------
        pcoords = np.array([[0.0], [0.9], [9.5]], np.float32)
        assign = bin_mapper.assign(pcoords, add_bins=True)
        assert list(assign) == [0, 0, 1]
        assert bin_mapper.nbins == 3

        assign = bin_mapper.assign(pcoords)
        assert list(assign) == [0, 0, 2]

        #---------
        pcoords = np.array([[1.2], [10.6]], np.float32)
        assign = bin_mapper.assign(pcoords, add_bins=True)
        assert list(assign) == [0, 2]
        assert bin_mapper.nbins == 5

        assign = bin_mapper.assign(pcoords)
        assert list(assign) == [3, 4]

