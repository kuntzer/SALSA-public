import numpy as np

def declination(x,y,z):
    return np.arctan2(z, np.sqrt(x*x+y*y) )

def right_ascension(x,y):
    return np.arctan2(y, x)

def R_3d(x,y,z):
    return np.sqrt(x*x+y*y+z*z)

def radec2cart(alpha,delta,r):
    x = r*np.cos(alpha)*np.cos(delta)
    y = r*np.sin(alpha)*np.cos(delta)
    z = r*np.sin(delta)
    return x,y,z

def radec2cart_vector(alpha,delta,r):
    x, y, z= (alpha,delta,r)
    return np.asarray([x,y,z])

def cart2radec(x,y,z):
    return right_ascension(x,y), declination(x,y,z), R_3d(x,y,z)

def rev(angle):
    twopi = 2.*np.pi
    angle = angle - int(angle/twopi)*twopi
    if (angle < 0):
        angle = angle + twopi
    return angle
''' Vectorisation of the function (for numpy.arrays) See http://docs.scipy.org/doc/numpy/reference/generated/numpy.vectorize.html '''
vrev = np.vectorize(rev)

def ROT_MATRIX_X(theta, x, y, z):
#     For the inverse matrix replace theta by -theta

    xx = x
    yy = y
    zz = z

    x = xx
    y = yy*np.cos(theta) - zz*np.sin(theta)
    z = yy*np.sin(theta) + zz*np.cos(theta)

    return x,y,z

def ROT_MATRIX_Y(theta, x, y, z):
#     For the inverse matrix replace theta by -theta
    xx = x
    yy = y
    zz = z

    x = xx*np.cos(theta) + zz*np.sin(theta)
    y = yy
    z = -xx*np.sin(theta) + zz*np.cos(theta)

    return x,y,z


def ROT_MATRIX_Z(theta, x, y, z):
#     For the inverse matrix replace theta by -theta
    xx = x
    yy = y
    zz = z
       
    x = xx*np.cos(theta) - yy*np.sin(theta)
    y = xx*np.sin(theta) + yy*np.cos(theta)
    z = zz

    return x,y,z

def limit_Earth(gamma,RA_SAT,DEC_SAT,n=100):
    def sign(x,y):
        if (x > 0 and y > 0) or (x < 0 and y < 0): return 1
        else: return -1

    LIMIT = np.zeros([n,2])

    x = np.cos(gamma)*np.cos(RA_SAT)
    y = np.cos(gamma)*np.sin(RA_SAT)
    z = sign(1.,DEC_SAT)*np.sin(gamma)

    theta = -RA_SAT
    x, y, z = ROT_MATRIX_Z(theta, x, y, z)

    theta = DEC_SAT
    x, y, z = ROT_MATRIX_Y(theta, x, y, z)

    xx = x
    yy = y
    zz = z

    for i in range(1, n+1):
#        print i, x, y, z
        x = xx
        y = yy
        z = zz

        theta = float(i-1)*2.*np.pi/float(n-1)

        x,y,z = ROT_MATRIX_X(theta, x, y, z) 
#        print i, x, y, z
        theta = -DEC_SAT
        x,y,z = ROT_MATRIX_Y(theta, x, y, z)
#        print i, x, y, z
        theta = RA_SAT
        x,y,z = ROT_MATRIX_Z(theta, x, y, z)
#        print i, x, y, z
        LIMIT[i-1,0]= rev(np.arctan2(y,x))
        LIMIT[i-1,1]= np.arctan2(z, np.sqrt(x*x+y*y))
    return LIMIT

def limit_theta2(RA_SAT,DEC_SAT,angle,n=300):
    theta = np.linspace(0,2*np.pi,n)
    LIMIT = np.zeros([n,2])

    LIMIT[:,0] = vrev((np.cos(theta))*angle+RA_SAT)
    LIMIT[:,1] = np.sin(theta)*angle+DEC_SAT

#    import pylab as plt
#    plt.figure()
#    plt.plot(LIMIT[:,0],LIMIT[:,1])
#    plt.show()

    return LIMIT

def SphericalDistance(ra1,dec1,ra2,dec2):
    '''Compute the distance on the great circle. Computations are based upon http://en.wikipedia.org/wiki/Great-circle_distance#Formulas'''
    from numpy import sin, cos, arctan2
    ls = ra1
    lf = ra2

    ps = dec1
    pf = dec2

    dl = lf-ls
    dp = pf-ps

    num1 = cos(pf)*sin(dl)
    num2 = (cos(ps)*sin(pf) - sin(ps)*cos(pf)*cos(dl))

    den1 = sin(ps)*sin(pf)
    den2 = cos(ps)*cos(pf)*cos(dl)

    num = np.sqrt(num1 * num1 + num2*num2)
    den = den1 + den2

    return arctan2(num,den)
vSphericalDistance = np.vectorize(SphericalDistance)

def rev_d(angle):
    if angle < -np.pi/2: angle += np.pi
    elif angle > np.pi/2: angle -= np.pi
    return angle
vrev_d = np.vectorize(rev_d)
