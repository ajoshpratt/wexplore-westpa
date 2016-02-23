# Courtesy of Karl Debiec
# See https://chong.chem.pitt.edu/wewiki/Introductory_Tutorial:_GROMACS
import h5py, numpy, pylab
from matplotlib import pyplot as plt
fluxanl              = h5py.File('fluxanl.h5')
west                 = h5py.File('west.h5')
first_binding        = fluxanl['target_flux']['target_0']['flux_evolution']
x = []
current_walkers = 0
for i in range(0, first_binding['expected'].shape[0]):
    current_walkers += west['summary']['n_particles'][i]
    x.append(current_walkers*.5)
fig, ax = plt.subplots(1)
ax.plot(x, [0.046164739614646955]*len(x), color='black', lw=2)
ax.plot(x, [0.023633008512743959]*len(x), color='black', lw=2)
ax.fill_between(x, [0.023633008512743959]*len(x), [0.046164739614646955]*len(x), facecolor='black', alpha=0.5)
ax.plot(x, first_binding['expected']*2, color='red', lw=2)
ax.plot(x, first_binding['ci_lbound']*2, color='red', lw=2)
ax.plot(x, first_binding['ci_ubound']*2, color='red', lw=2)
ax.fill_between(x, first_binding['ci_lbound']*2, first_binding['ci_ubound']*2, facecolor='red', alpha=0.5)

plt.xlabel("Total Time")
plt.ylabel("Instantaneous Flux $(\\frac{1}{ps})$")
plt.axis([0,9000, 0, 0.05])
plt.show()

