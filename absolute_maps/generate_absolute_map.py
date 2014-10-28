import resources
import sys
import numpy as np
import os
import multiprocessing

lock=multiprocessing.Lock()

class AbsoluteObservationMap(object):
    
    def __init__(self, orbit_id, grid={"type":"cartesian","resolution":[80,40]}, outdir=None, verbose=True,ncpu=1):
        self.orbit_id = orbit_id
        print 'Absolute Observation map for orbit id %s' % self.orbit_id
        self.grid=grid
        self.generate_grid_points()
        self.verbose = verbose
        if outdir==None:
            self.outdir = "absolute_map_%s" % self.orbit_id
        else:
            self.outdir=outdir
            
        self.orbital_period=None
        self.do_monitor_SAA=False
        if ncpu==0:
            self.ncpu=multiprocessing.cpu_count()
        else: 
            self.ncpu=ncpu
        
    def generate_grid_points(self):
        """ Generates the grid points necessary for the computation """
        if self.grid["type"]=="cartesian":
            ra_i = 0.
            ra_f = 2*np.pi
            ra_step=(ra_f-ra_i)/self.grid["resolution"][0]
            dec_i = -np.pi/2.
            dec_f = np.pi/2.
            dec_step=(dec_f-dec_i)/self.grid["resolution"][1]
            
            ra = np.arange(ra_i+ra_step/2, ra_f, ra_step)
            dec= np.arange(dec_i+dec_step/2, dec_f, dec_step)
            self.grid["points"]=np.array([[a,d] for a in ra for d in dec]) 
            self.grid["mesh"] = np.meshgrid(ra,dec)
        else:
            raise ValueError("Grid type unknown")
        
    def set_constraints(self,angles={"sun":120,"moon":5,"earth":35},altitudes={"atmosphere":100e3},shutdowns=False):
        """ Set the observability constraints. Angles in degrees, altitudes in km. """
        [setattr(self, "angle_%s"%angle, np.deg2rad(angles[angle])) for angle in angles]
        [setattr(self, "altitude_%s"%altitude, altitudes[altitude]) for altitude in altitudes]
        
    def set_orbital_period(self,period=None, apogee=None, perigee=None):
        if not period==None and (apogee==None and perigee==None):
            self.orbital_period=period
        elif period==None and not (apogee==None and perigee==None):
            a=(apogee+perigee)/2.+resources.constants.R_earth
            self.orbital_period=2.*np.pi*np.sqrt(a*a*a/resources.constants.mu_Earth)
        else:
            raise ValueError("Orbital elements not understood.")
        
    def load_object_position(self,objects=["sun","earth","moon"],coord="latlon"):
        """ Loads every physical objects needed for the computation and saves what was loaded """
        if self.verbose:
            print "Loading positions of objects: ", objects
            sys.stdout.write("This can take a while... ")
            sys.stdout.flush() 
        for body in objects:
            angles, cart=resources.utils.load_position_file(self.orbit_id, body,coord=coord)
            setattr(self, "position_%s"%body, angles)
            setattr(self, "position_cart_%s"%body, cart)
        self.loaded_objects=objects
        if self.verbose: print "Done."
        
    def info(self):
        """
        Prints out all the variables for a given absolute map of the class
        """
        import inspect
    
        message = "All variables available for absolute map %s" % self.orbit_id        
        print message
        print '-'*len(message)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if (a[0].startswith('__') and a[0].endswith('__')): continue
            print a[0], "=", a[1]
            
    def test_angle(self,points,body,minute):
        """ Test points according to the rule of a given body """
        angular_diameters={"sun":0.545/180.*np.pi,"moon":0.5685/180.*np.pi} # must be in rad
        
        
        coord_body = getattr(self,"position_%s"%body)[minute]
        if body in angular_diameters: corr = angular_diameters[body]
        else: corr = 0.
        dist = np.abs(resources.maps.vSphericalDistance(points[:,0],points[:,1],coord_body[1],coord_body[2])-corr)
        dist[dist<=getattr(self,"angle_%s"%body)]=False
        dist[dist>getattr(self,"angle_%s"%body)]=True

        return np.asarray(dist,dtype=bool)
    
    def monitor_SAA(self,radiation_limit=2):
        pos=np.loadtxt(os.path.join('orbits',self.orbit_id,'spenvis_pos.txt'),delimiter=",",comments="'")
        flux=np.loadtxt(os.path.join('orbits',self.orbit_id,'spenvis_tpo.txt'), delimiter=',',comments='\'')
        assert np.shape(pos)[0]==np.shape(flux)[0]
        
        self.shutdowns_SAA = np.asarray([pos[:,2]/180.*np.pi,pos[:,1]/180.*np.pi,flux[:,-1]]).T
        self.do_monitor_SAA=True
        self.radiation_limit=radiation_limit
        self.shutdowns_SAA[self.shutdowns_SAA[:,2]<-1,2]=0
        self.shutdowns_SAA[self.shutdowns_SAA[:,2]<=self.radiation_limit,2]=False
        self.shutdowns_SAA[self.shutdowns_SAA[:,2]>self.radiation_limit,2]=True
            
    def compute_observability(self,start,end,time_step=1,filename='orbit_',show=False):
        """ Does the heavy lifting, computes for every time step the absolute observability map. 
        Can also skip orbit_step to speed the computations up. """
        if end==-1: end = np.shape(self.position_earth)[0]
        
        self.filename = filename
        self.show = show
        self.time_step=time_step
        
        if not os.path.exists(self.outdir): os.mkdir(self.outdir)
        else: 
            print 'Stop. The folder already exists.'; exit();
        
        constraints=open(os.path.join(self.outdir,'constraints_%s.dat' % self.orbit_id),'w')
        print >> constraints, 'Exclusions angles\n==============='
        for body in self.loaded_objects:
            print >> constraints, body,'=',np.rad2deg(getattr(self,"angle_%s"%body)), 'deg'
            
        print >> constraints, '\nAltitudes\n==============='
        print >> constraints, 'Atmosphere=',self.altitude_atmosphere, 'm'
        constraints.close()
            
        if self.do_monitor_SAA: 
            file_SAA = open(os.path.join(self.outdir,'SAA_table_%s.dat' % (self.orbit_id)),'w')
            pos_sat = np.loadtxt(os.path.join('orbits',self.orbit_id,'sat.dat'))
            in_saa=False
        
        #latest_orbit_saved=0
        minute_table = open(os.path.join(self.outdir,'minute_table_%s.dat' % (self.orbit_id)),'w')
        minute_start=start

        latest_orbit_saved=0
        for minute in range(start,end,time_step):
            if float(minute)>(latest_orbit_saved+1)*self.orbital_period/60.:

                print >> minute_table, '%d,%d,%d' % (latest_orbit_saved+1,0,minute_start)
                print >> minute_table, '%d,%d,%d' % (latest_orbit_saved+1,minute-minute_start,minute-1)

                minute_start=minute
                latest_orbit_saved+=1
                
                sys.stdout.write('\rPreparing orbit %d' % (latest_orbit_saved))
                sys.stdout.flush()
                
            if self.do_monitor_SAA:
                lon_id=resources.utils.find_nearest(self.shutdowns_SAA[:,1], pos_sat[minute,2])
                #print np.round(pos_sat[minute,1]*180/np.pi), np.round(pos_sat[minute,2]*180/np.pi)
                lon = self.shutdowns_SAA[lon_id,1]
                try :
                    _=np.shape(lon)[0]
                    lon=lon[0]
                except: pass
                lat = self.shutdowns_SAA[lon==self.shutdowns_SAA[:,1]]
                is_saa_active=resources.utils.find_nearest(lat[:,0], pos_sat[minute,1])
                #print lat[is_saa_active]*180/np.pi
                if not in_saa==lat[is_saa_active,2]:
                    #lock.acquire()
                    print >> file_SAA, '%d,%d' % (minute, lat[is_saa_active,2]) 
                    #lock.release()
                in_saa=lat[is_saa_active,2]
                
        minute_table.close()
        if self.do_monitor_SAA: file_SAA.close()
                
        print 'Starting computations on %d cpu' % self.ncpu
        
        from itertools import repeat
        mpo=self.orbital_period/60.
        orbit_end = np.floor(end/mpo)
        orbit_start = np.floor(start/mpo)
        orbits_per_cpu=(orbit_end-orbit_start)/self.ncpu
        orbits_per_cpu=np.floor(orbits_per_cpu).astype(int)

        starts=[]
        ends=[]
        orbit_ini=[]
        for n in range(self.ncpu):
            if n==0: start=start
            else: start = n*orbits_per_cpu*mpo+1
            end = (n+1)*mpo*orbits_per_cpu
            starts.append(int(start))
            ends.append(int(end))
            orbit_ini.append(n*orbits_per_cpu)

        ends[-1]=int(np.ceil(end))

        inputs=zip(repeat(self),starts,ends,orbit_ini)
        pool = multiprocessing.Pool(processes=self.ncpu)
        #inputs=
        pool.map(_worker, inputs)
        pool.close()
        pool.join()
            
    def show(self, proj='moll', lon_0=180, tmap=None, coord=None):
        from mpl_toolkits.basemap import Basemap
        import pylab as plt
        import resources.figures as figures
        figures.set_fancy()
        if coord==None:
            ra = np.rad2deg(self.grid['points'][:,0])
            dec = np.rad2deg(self.grid['points'][:,1])
        else:
            ra = np.rad2deg(coord[:,0])
            dec = np.rad2deg(coord[:,1])
        
        fig = plt.figure()
        m = Basemap(projection=proj,lon_0=lon_0)#,celestial=True) # celestial=True inverses alpha (East towards the right)
        m.drawparallels(np.arange(-60.,90.,30.),labels=[1,0,0,0])
        m.drawmeridians(np.arange(0.,360.,30.))
        
        ra__ = np.arange(0., 360., 30.)
        x, y = m(ra__,ra__*0)
        for x,y,t in zip(x,y,ra__):
            plt.text(x, y, figures.format_degree(t), color='black', ha='center', weight='black', size='small') ##93c6ed
        if tmap==None:
            m.scatter(ra,dec,latlon=True,marker='x',s=20,cmap=plt.cm.binary)
        else:
            m.scatter(ra,dec,c=tmap,latlon=True,marker='x',s=20,cmap=plt.cm.binary)
        plt.show()
        
def _worker(inputs):
    """
    Worker function that processes one _AbsoluteObservationMap :compute_observability: object.
    """
    obj, start, end, latest_orbit_saved = inputs

    obs_coord=None
    sys.stdout.write('\rComputing orbit %d' % (latest_orbit_saved))
    sys.stdout.flush()
    for minute in range(start,end+obj.time_step,obj.time_step):
        if float(minute)>(latest_orbit_saved+1)*obj.orbital_period/60.:
            np.savetxt(os.path.join(obj.outdir,'%s%d.dat' % (obj.filename,latest_orbit_saved+1)),obs_coord,fmt='%3.6f %3.6f %3.6f')
            obs_coord=None
            minute_start=minute
            latest_orbit_saved+=1
            
            sys.stdout.write('\rComputing orbit %d' % (latest_orbit_saved))
            sys.stdout.flush()
        

        # Alternate method. Slower as we compute every thing again!
        #for ii, body in enumerate(obj.loaded_objects):
        #    if ii==0: 
        #        tmap=np.ones_like(obj.grid['points'][:,0],dtype=bool)
        #    if body=="earth":
        #        ...
        #    else:
        #        tmap = np.bitwise_and(tmap, obj.test_angle(obj.grid["points"],body,minute))
        #############################
        coord=obj.grid['points']
        
        for ii, body in enumerate(obj.loaded_objects):
            if ii==0: 
                tmap=np.ones_like(coord[:,0],dtype=bool)
            
            if body=="earth":
                # 1. Are we eclipsed by the Earth ?
                _,ra_e,dec_e,r_e = getattr(obj,"position_%s"%body)[minute]
                try: atmosphere=obj.altitude_atmosphere
                except: atmosphere=0.
                # This is the angle at which we see the Earth.
                alpha=np.arcsin((resources.constants.R_earth+atmosphere)/r_e/1e3) 
                # We test the angle to see whether we are not observing "through" the Earth
                tmap = resources.maps.vSphericalDistance(coord[:,0],coord[:,1],ra_e, dec_e)>=alpha
                # and we remove those which are eclipsed
                coord=coord[tmap,:]
                
                # reinitialise the truth map
                #coord=np.array([[0,0],[0,0]])
                tmap=np.ones_like(coord[:,0],dtype=bool)
                # 2. Check the SL exclusion angle (i.e. angle_earth) 
                # This technique is flawed... But I couldn't find a better and still fast solution (otherwise have to compute as in SL code...) ==> MERGE TWO CODES??
                # See http://arxiv.org/abs/1310.7800, appendix 1
                # This 2. part is directly taken from the code of Luzius Kronig.
                # calculate gracing point (without atmosphere)
                sun_pos=obj.position_cart_sun[minute]
                sun_pos=sun_pos[1:]
                #earth_pos=resources.maps.radec2cart_vector(ra_e,dec_e,r_e)
                earth_pos = getattr(obj,"position_cart_%s"%body)[minute]
                earth_pos = earth_pos[1:]
                earth_dir=earth_pos/r_e
                
                alpha2=np.arccos(resources.constants.R_earth/r_e/1e3)

                star_pos=np.zeros([np.shape(coord)[0],3])
                star_pos[:,0]=np.cos(coord[:,1])*np.cos(coord[:,0])
                star_pos[:,1]=np.cos(coord[:,1])*np.sin(coord[:,0])
                star_pos[:,2]=np.sin(coord[:,1])
                b=np.cross(earth_dir,star_pos)
                
                a=np.cross(b,earth_dir)
                a = np.asarray([a[ii,:]/nn for ii, nn in enumerate(resources.utils.norm(a))]) 
                gracing_pos2=earth_pos-resources.constants.R_earth/1.e3*np.cos(alpha2)*earth_dir+resources.constants.R_earth/1.e3*np.sin(alpha2)*a

                # okay until here
                
                # calculate terminator point and terminator vector from satellite
                b=np.cross((sun_pos-earth_pos),b)
                b=np.asarray([b[ii,:]/nn for ii, nn in enumerate(resources.utils.norm(b))])
                
                terminator_pos=earth_pos+resources.constants.R_earth*b/1.e3
                """print
                print 'earth pos', earth_pos
                print earth_pos-resources.constants.R_earth/1.e3*np.cos(alpha2)*earth_dir
                print 'star pos', star_pos
                print gracing_pos2
                exit()"""
                
                earth_sun_vector=sun_pos-earth_pos
                earth_gracing_vector=gracing_pos2-earth_pos
                
                #norm_egracing = np.sqrt((earth_gracing_vector*earth_gracing_vector).sum(axis=1))
                earth_gracing_vector= np.asarray([earth_gracing_vector[ii,:]/nn for ii, nn in enumerate(resources.utils.norm(earth_gracing_vector))])
                #UNUSED:dot_product = [np.dot(earth_sun_vector,v) for v in earth_gracing_vector]
                norm_earth_sun_vector = np.sqrt((earth_sun_vector*earth_sun_vector).sum())
                
                #print earth_sun_vector
                #print earth_gracing_vector[5,]
                #print np.dot(earth_sun_vector,earth_gracing_vector[5]), '<<< dot'
                #print norm_earth_sun_vector
                #test_angle=acos(dot(earth_sun_vector,earth_gracing_vector)/(norm(earth_sun_vector)*norm(earth_gracing_vector)));
                 
                #test_angles=np.arccos(dot_product /(norm_earth_sun_vector)) ### REDO THIS ! It STILL DOESNT WORK
                #print test_angles[5]
                #exit()
                                    
                # TODO: The following can be improved to a compute the angles for the gracing_pos all at once and comparing the numpy array => faster.
                for i, [c_a,c_d] in enumerate(coord):
                    test_angle=np.arccos(np.dot(earth_sun_vector,earth_gracing_vector[i,:])/(norm_earth_sun_vector))
                        
                    if test_angle<np.pi/2:#[i] < np.pi/2.:
                        
                        ra_t, dec_t,_ = resources.maps.cart2radec(gracing_pos2[i,0], gracing_pos2[i,1], gracing_pos2[i,2])
                        #print minute, test_angle[i]*180./np.pi,resources.maps.vSphericalDistance(c_a,c_d,ra_t, dec_t)*180/np.pi,
                        #ra_t, dec_t,_ = resources.maps.cart2radec(terminator_pos[i,0], terminator_pos[i,1], terminator_pos[i,2])
                        #print resources.maps.vSphericalDistance(c_a,c_d,ra_t, dec_t)*180/np.pi
                    else: 
                        ra_t, dec_t,_ = resources.maps.cart2radec(terminator_pos[i,0], terminator_pos[i,1], terminator_pos[i,2])
                    straylight_angle=resources.maps.vSphericalDistance(c_a,c_d,ra_t, dec_t)
                    
                    if straylight_angle>=getattr(obj,"angle_%s"%body): tmap[i]=True
                    else: tmap[i]=False
                    
                # and we remove those which are > SL angle
                coord=coord[tmap,:]
                
            else:
                # For the other bodies, we simply take a fixed exclusion angle.
                tmap = obj.test_angle(coord,body,minute)
                coord=coord[tmap,:]
                
        if obj.show:obj.show(tmap=None,coord=coord)
            
        obs_coord_tmp=np.array([np.ones_like(coord[:,0],dtype=int)*minute, coord[:,0], coord[:,1]]).T
        if obs_coord==None:
            obs_coord=obs_coord_tmp
        else:
            obs_coord=np.vstack([obs_coord,obs_coord_tmp])
        del coord
        
    np.savetxt(os.path.join(obj.outdir,'%s%d.dat' % (obj.filename,latest_orbit_saved+1))
               ,obs_coord,fmt='%3.6f %3.6f %3.6f')
    del obs_coord

if __name__=="__main__":
    # Put the same name here as for the folder
    altitude=800
    cheops=AbsoluteObservationMap("SSO%d" % altitude,
                                  grid={"type":"cartesian","resolution":[80,40]},
                                  outdir="absolute_map_800_80x40",ncpu=0)
    # This is useful for saving the files.
    cheops.set_orbital_period(apogee=altitude*1e3, perigee=altitude*1e3)
    # We want to monitor if is the SAA is crossed:
    cheops.monitor_SAA(radiation_limit=2)
    # Exclusion angles in degrees are set. For the earth, you must specify a value. 
    # Altitude of the atmosphere in m.
    cheops.set_constraints(angles={"sun":120,"moon":5,"earth":35},
                           altitudes={"atmosphere":100e3})
    # Loading the corresponding file. The order here matters as it will be the order the potential 
    # targets are removed in. I recommend Sun, Moon and then Earth
    # as it is the most time consuming part. objects=["sun","moon","earth"]
    cheops.load_object_position(objects=["sun","moon","earth"],coord="cartesian")
    # Specify when to start and when to finish in unit of the time resolution in the
    # position file [For CHEOPS==1 minute]. 
    # if end=-1 compute for every available timestep
    cheops.compute_observability(start=0, end=-1)
    # Prints everything in the object cheops.
    #print absmap.info()