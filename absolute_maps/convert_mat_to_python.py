import scipy.io
import numpy as np

alt=620
mat=scipy.io.loadmat('orbits/SSO%d.mat' % alt)

print 'List of all variable in file:'
for m in mat:
    print m
print '--'*20

lat=mat['lat'][0]
t=mat['time']
lon=mat['lon'][0]
print 'Converting geographical pos of sat'
f = open('orbits/SSO%d/sat.dat'% alt,'w')
print >> f, '# t, lon, lat'
for a,b,c in zip(t,lon,lat):
    print >> f, a[0],b,c
f.close()
#np.savetxt('orbits/SSO700/sat_latlon.dat',tosave)

del lat,lon

print 'Converting pos of sun'
sun=mat['sun_pos']
f = open('orbits/SSO%d/sun.dat'% alt,'w')
print >> f, '# t [min], x,y,z [km]'
for ti,[a,b,c] in zip(t,sun):
    print >> f, ti[0], a,b,c
f.close()

del sun

print 'Converting pos of Earth'
earth_pos=mat['earth_pos']
f = open('orbits/SSO%d/earth.dat'% alt,'w')
print >> f, '# t [min], x,y,z [km]'
for ti,[a,b,c] in zip(t,earth_pos):
    print >> f, ti[0], a,b,c
f.close()

print 'Converting pos of orbit'
f = open('orbits/SSO%d/orbit.dat'% alt,'w')
print >> f, '# t [min], x,y,z [km]'
for ti,[a,b,c] in zip(t,earth_pos):
    print >> f, ti[0], -a,-b,-c
f.close()

del earth_pos
print 'Converting pos of moon'
moon=mat['moon_pos']
f = open('orbits/SSO%d/moon.dat'% alt,'w')
print >> f, '# t [min], x,y,z [km]'
for ti,[a,b,c] in zip(t,moon):
    print >> f, ti[0], a,b,c
f.close()


print 'Conversion done.'