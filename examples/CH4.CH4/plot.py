# Courtesy of Karl Debiec
# See https://chong.chem.pitt.edu/wewiki/Introductory_Tutorial:_GROMACS
import h5py, numpy, pylab
fluxanl              = h5py.File('fluxanl.h5')
flux                 = numpy.zeros(100)
first_binding        = fluxanl['target_flux']['target_0']['flux_evolution']
pylab.plot(first_binding['expected']*2, color='black')
pylab.plot(first_binding['ci_lbound']*2, color='grey')
pylab.plot(first_binding['ci_ubound']*2, color='grey')
pylab.xlabel("Iteration")
pylab.ylabel("Instantaneous Flux $(\\frac{1}{ps})$")
pylab.show()

