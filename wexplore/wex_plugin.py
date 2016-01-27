from __future__ import division; __metaclass__ = type
import logging
log = logging.getLogger(__name__)

import numpy as np
import time

import westpa, west
from westpa import extloader
from westpa.yamlcfg import check_bool, ConfigItemMissing
from westpa.binning import RectilinearBinMapper, RecursiveBinMapper, FuncBinMapper
from wexplore import wexplore
import cStringIO
from west.states import TargetState, InitialState, BasisState
from west import Segment

class WExploreDriver(object):
    def __init__(self, sim_manager, plugin_config):
        super(WExploreDriver, self).__init__()

        if not sim_manager.work_manager.is_master:
                return

        self.sim_manager = sim_manager
        self.data_manager = sim_manager.data_manager
        self.system = sim_manager.system
        self.we_driver = sim_manager.we_driver
        self.priority = plugin_config.get('priority', 0)

        self.init_from_data = check_bool(plugin_config.get('init_from_data', True))

        self.max_replicas = self.system.max_replicas

        # Initialize bin mapper from data in h5 file if available
        bin_mapper = self.init_bin_mapper()

        print(bin_mapper)
        if bin_mapper:
            self.system.bin_mapper = bin_mapper
            self.we_driver.bin_mapper = bin_mapper

        # Register callback
        sim_manager.register_callback(sim_manager.pre_we, self.pre_we, self.priority)
        sim_manager.register_callback(sim_manager.run_we, self.run_we, self.priority)
        sim_manager.register_callback(sim_manager.pre_prepare_iteration, self.target_counts, self.priority)

    def init_bin_mapper(self):
        self.data_manager.open_backing()

        with self.data_manager.lock:
            #n_iter = max(self.data_manager.current_iteration - 1, 1)
            n_iter = max(self.data_manager.current_iteration, 1)
            iter_group = self.data_manager.get_iter_group(n_iter)

            # First attempt to initialize binmapper from data rather than system
            bin_mapper = None
            if self.init_from_data:
                log.info('Attempting to initialize WEX bin mapper from data')

                try:
                    binhash = iter_group.attrs['binhash']
                    bin_mapper = self.data_manager.get_bin_mapper(binhash)

                except:
                    log.warning('Initializing bins from data failed; Using definition in system instead.')
                    centers = self.system.bin_mapper.centers
            else:
                log.info('Initializing bin mapper from system definition')

        self.data_manager.close_backing()

        return bin_mapper

    def pre_we(self):
        starttime = time.time()
        segments = self.sim_manager.segments.values()
        bin_mapper = self.system.bin_mapper
        wexploreMappers = []
        for i in self.system.bin_mapper.mapper_list:
            if isinstance(i['base_mapper'], wexplore.WExploreBinMapper) == True:
                wexploreMappers.append(i)
        wexploreMapper = wexploreMappers[0]['base_mapper']

        final_pcoords = np.empty((len(segments), self.system.pcoord_ndim), dtype=self.system.pcoord_dtype)

        for iseg, segment in enumerate(segments):
            final_pcoords[iseg] = segment.pcoord[-1,:]

        hash_init = wexploreMapper.last_graph
        kwargs = {'add_bins': True}
        assignments = bin_mapper.assign(final_pcoords, **kwargs)

        # Re-assign segments to new bins if bin_mapper changes
        if wexploreMapper.last_graph != hash_init:
            # Redo target states.  Actually, this may need to happen later?
            target_states = []
            #self.system.tstates = ['unbound,5.0,1.0']
            tstates = self.system.tstates
            tstates_strio = cStringIO.StringIO('\n'.join(tstates).replace(',', ' '))
            target_states.extend(TargetState.states_from_file(tstates_strio, self.system.pcoord_dtype))
            self.data_manager.save_target_states(target_states)
            self.sim_manager.report_target_states(target_states)

            # I think we need to force the mapper to reupdate...
            bin_mapper = self.system.initialize_mappers(wexploreMapper)
            # Now, update the local bin mapper.
            self.we_driver.bin_mapper = bin_mapper
            self.sim_manager.bin_mapper = bin_mapper


            # Reset we_driver internal data structures
            initial_states = []
            for key,value in self.we_driver.avail_initial_states.iteritems():
                initial_states.append(value)
            self.we_driver.new_iteration(target_states=target_states, initial_states=initial_states)

            # Re-assign segments
            self.we_driver.assign(segments)

        # Get assignments. Should use cached assignments
        # We'll pull this from the system, now...
        #start_index = self.system.bin_mapper.mapper_list[1]['start_index']
        start_index = wexploreMappers[0]['mapper'].start_index
        print("This is the start index! " + str(start_index))
        self.target_counts(coords=final_pcoords)

        endtime = time.time()

        # Report level statistics
        s = 1
        westpa.rc.pstatus('--wexplore-stats--------------------')
        westpa.rc.pstatus('wallclock time: {:.3f} s'.format(endtime - starttime))
        westpa.rc.pstatus('')
        for li, level in enumerate(wexploreMapper.level_indices):
            s *= wexploreMapper.n_regions[li]
            westpa.rc.pstatus('Level {}: {} cells ({} max)'.format(li, len(level), s))
        westpa.rc.pstatus('------------------------------------')
        westpa.rc.pflush()

    def run_we(self):
        '''Run the weighted ensemble algorithm based on the binning in self.final_bins and
        the recycled particles in self.to_recycle, creating and committing the next iteration's
        segments to storage as well.'''
        
        # The WE driver now does almost everything; we just have to record the
        # mapper used for binning this iteration, and update initial states
        # that have been used
         
        try:
            pickled, hashed = self.we_driver.bin_mapper.pickle_and_hash()
        except PickleError:
            pickled = hashed = ''

        self.bin_mapper_hash = hashed

        self.construct_next()
        
        if self.we_driver.used_initial_states:
            for initial_state in self.we_driver.used_initial_states.itervalues():
                initial_state.iter_used = self.sim_manager.n_iter+1
            self.data_manager.update_initial_states(self.we_driver.used_initial_states.values())
            
        self.data_manager.update_segments(self.sim_manager.n_iter,self.sim_manager.segments.values())
        
        self.data_manager.require_iter_group(self.sim_manager.n_iter+1)
        self.data_manager.save_iter_binning(self.sim_manager.n_iter+1, hashed, pickled, self.we_driver.bin_target_counts)
        
        # Report on recycling
        recycling_events = {}
        for nw in self.we_driver.new_weights:
            try:
                recycling_events[nw.target_state_id].add(nw.weight)
            except KeyError:
                recycling_events[nw.target_state_id] = set([nw.weight])
        
        tstates_by_id = {state.state_id: state for state in self.we_driver.target_states.itervalues()}
        for tstate_id, weights in recycling_events.iteritems():
            tstate = tstates_by_id[tstate_id]
            self.sim_manager.rc.pstatus('Recycled {:g} probability ({:d} walkers) from target state {!r}'.format(sum(weights),
                                                                                                     len(weights),
                                                                                                     tstate.label))

    def construct_next(self):
        '''Construct walkers for the next iteration, by running weighted ensemble recycling
        and bin/split/merge on the segments previously assigned to bins using ``assign``.
        Enough unused initial states must be present in ``self.avail_initial_states`` for every recycled
        walker to be assigned an initial state.
        
        After this function completes, ``self.flux_matrix`` contains a valid flux matrix for this
        iteration (including any contributions from recycling from the previous iteration), and
        ``self.next_iter_segments`` contains a list of segments ready for the next iteration,
        with appropriate values set for weight, endpoint type, parent walkers, and so on.        
        '''
        
        self.we_driver._prep_we()
        
        # Create new segments for the next iteration        
        # We assume that everything is going to continue without being touched by recycling or WE, and
        # adjust later
        new_pcoord_array = self.system.new_pcoord_array
        n_iter = None
        segments = []
                
        for ibin, _bin in enumerate(self.we_driver.final_binning):
            for segment in _bin:
                if n_iter is None:
                    n_iter = segment.n_iter
                else:
                    assert segment.n_iter == n_iter
                    
                segment.endpoint_type = Segment.SEG_ENDPOINT_CONTINUES
                new_segment = Segment(n_iter=segment.n_iter+1,
                                      parent_id=segment.seg_id,
                                      weight=segment.weight,
                                      wtg_parent_ids=[segment.seg_id],
                                      pcoord=new_pcoord_array(),
                                      status=Segment.SEG_STATUS_PREPARED)
                new_segment.pcoord[0] = segment.pcoord[-1]
                self.we_driver.next_iter_binning[ibin].add(new_segment)

                
                # Store a link to the parent segment, so we can update its endpoint status as we need,
                # based on its ID
                self.we_driver._parent_map[segment.seg_id] = segment
                
        self.we_driver._recycle_walkers()
        wexploreMappers = []
        for i in self.system.bin_mapper.mapper_list:
            if isinstance(i['base_mapper'], wexplore.WExploreBinMapper) == True:
                wexploreMappers.append(i)
        wexploreMapper = wexploreMappers[0]['base_mapper']
        self.system.initialize_mappers(wexploreMapper)
        self.target_counts()
        self.we_driver._run_we()
        
        log.debug('used initial states: {!r}'.format(self.we_driver.used_initial_states))
        log.debug('available initial states: {!r}'.format(self.we_driver.avail_initial_states))

    def target_counts(self, coords=None):
        inmapper = []
        segments = []
        #for i in self.we_driver.next_iter_binning:
        if coords == None:
            try:
                for i in self.we_driver.next_iter_segments:
                    segments.append(i)
                final_pcoords = np.empty((len(segments), self.system.pcoord_ndim), dtype=self.system.pcoord_dtype)

                for iseg, segment in enumerate(segments):
                    final_pcoords[iseg] = segment.pcoord[0]
            except(RuntimeError):
                # This implies we're calling it to restart an iteration.  We just need to set the appropriate list length.
                # This is to handle any target states present; when a simulation is restarted, it resets the target_bin
                # count to 0.  If it initializes this from the system, the list is generally not long enough.
                # Setting the counts to 1 shouldn't affect the running simulation.
                target_counts = [1] * self.system.bin_mapper.nbins
                self.system.bin_target_counts = target_counts
                self.we_driver.bin_target_counts = target_counts
                return 0

        else:
            final_pcoords = coords

        assignments = self.system.bin_mapper.assign(final_pcoords)
        wexploreMappers = []
        for i in self.system.bin_mapper.mapper_list:
            if isinstance(i['base_mapper'], wexplore.WExploreBinMapper) == True:
                wexploreMappers.append(i)
        wexploreMapper = wexploreMappers[0]['base_mapper']
        start_index = wexploreMappers[0]['mapper'].start_index
        # Only assign, to the wexploreMapper, those which should fall into it normally.
        for ii,i in enumerate(assignments):
            if i >= start_index and i < start_index + wexploreMapper.nbins:
                inmapper.append(final_pcoords[ii])
        assignmentswexplore = wexploreMapper.assign(inmapper)
        target_counts = wexploreMapper.balance_replicas(self.max_replicas, assignmentswexplore)
        total_bins = 0
        for imapper,mapper in enumerate(self.system.bin_mapper.mapper_list):
            total_bins += mapper['mapper'].nbins
        old_list = [0] * total_bins
        for imapper,mapper in enumerate(self.system.bin_mapper.mapper_list):
            start = mapper['start_index']
            nbins = mapper['base_mapper'].nbins
            targetc = mapper['base_mapper'].bin_target_counts
            old_list[start:start+nbins] = [targetc] * nbins
            
        old_list[start_index:start_index+wexploreMapper.nbins] = target_counts
        target_counts = old_list

        self.system.bin_target_counts = target_counts
        self.we_driver.bin_target_counts = target_counts
