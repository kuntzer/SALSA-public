import resources.figures as figures
import numpy as np
import pylab as plt
# WARNING: This is an obselete version of the code...

pst = np.loadtxt('straylight_orbitID_p/INPUT/pst.dat')
print pst
figures.set_fancy()
fig=plt.figure()
plt.xlabel(r'$\theta$')
plt.ylabel('PST')
plt.semilogy(pst[:,0],pst[:,1], lw=2)
plt.plot([35,35],[np.amin(pst[:,1]),np.amax(pst[:,1])],lw=3,color='r')
plt.grid()


figures.savefig('orbitID_figures/pst',fig,fancy=True)
plt.show()
