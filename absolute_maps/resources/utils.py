import numpy as np
import os

def load_position_file(ID,body,coord,path='orbits/'):
    path=os.path.join(path,ID,'%s.dat' % body)
    pos = np.loadtxt(path)
    if coord=='latlon':
        return pos
    elif coord=='cartesian':
        import resources.maps as maps
        alpha = maps.right_ascension(pos[:,1],pos[:,2])
        delta = maps.declination(pos[:,1],pos[:,2],pos[:,3])
        r = maps.R_3d(pos[:,1],pos[:,2],pos[:,3])
        return np.array([pos[:,0],alpha,delta,r]).T, pos
    else:
        raise NotImplemented()
    
def norm(b):
    return np.sqrt((b*b).sum(axis=1))

def find_nearest(array,value):
    ''' Find nearest value is an array '''
    idx = (np.abs(array-value)).argmin()
    return idx