c     ==========================================================================

	double precision function SphericalDistance(RA_SUN,DEC_SUN,RA,DEC)
c 	Compute the distance on the great circle. 
c 	Computations are based upon 
c	http://en.wikipedia.org/wiki/Great-circle_distance#Formulas
c	This is the *Angular* spherical distance in radians

	implicit none
	double precision RA_SUN, DEC_SUN, RA, DEC
	double precision ls, lf, ps, pf, dl, dp, num1, num2, den1, den2
	double precision num, den, angular_distance

	ls = RA_SUN
	lf = RA

	ps = DEC_SUN
	pf = DEC

	dl = lf-ls
	dp = pf-ps

	num1 = dcos(pf) * dsin(dl)
	num2 = dcos(ps) * dsin(pf) - dsin(ps)*dcos(pf)*dcos(dl)

	den1 = dsin(ps) * dsin(pf)
	den2 = dcos(ps) * dcos(pf) * dcos(dl)

	num = dsqrt(num1 * num1 + num2 * num2)
	den = den1 + den2

	SphericalDistance = datan2(num,den)

	return

	end

c     ==========================================================================
      subroutine ILLUMINATED_EARTH(RA_SUN, DEC_SUN, !Input
     &                             day_night)       !Output
c     Calculate the equatorial coordinates of the illuminated Earth
c	Compute which part of the Earth is illuminated and with what intensity 
c	(in arbitrary unit between 0 (total night) and 1 (zenith)) and outputs
c	the results in a table of the size of Earth's grid (n_alpha x n_delta).

      implicit none
      include 'constants'
      include 'parameters'

	double precision RA_SUN, DEC_SUN
	double precision day_night(1:n_alpha, 1:n_delta)
	double precision SphericalDistance
	double precision alpha_grid, delta_grid, dist
      integer k, j

	do k = 1, n_alpha
	  alpha_grid = 2.d0*pi*dble(k-1)/dble(n_alpha-1)
	  do j = 1, n_delta
	    delta_grid = -pi/2.d0+dble(j-1)*pi/dble(n_delta-1)

c	    This is the *Angular* spherical distance
	    dist = SphericalDistance(RA_SUN,DEC_SUN,alpha_grid,delta_grid)

	    if (dist <= pi/2.) then
c		Compute correct illumination of a given cell on the Earth's surface
		if (gradient) then
		  day_night(k,j) = day * dcos(dist)
		else
		  day_night(k,j) = day
		endif
	    else
c		Initialise the rest of the table.
		day_night(k,j) = night
	    endif

	  enddo
	enddo

	if (log_all_data) then
	  call write_illuminated_earth(day_night)
	endif

      return
      end
      
c     ========================================================================== 
      subroutine EARTH_SEEN_FROM_SAT (SAT, RA_SAT, DEC_SAT, day_night,	!Input
     &                                coord_illum) 
c	Returns the coordinates of the cells of the Earths that are 1) in daylight
c	and 2) visible the the satellite.               	!Output

      implicit none
      include 'constants'
      include 'parameters'

      double precision SAT(1:3), RA_SAT, DEC_SAT, R_SAT 
	double precision alpha_grid, delta_grid
      double precision alpha, SphericalDistance, dist
      double precision RA_LIMIT(1:n_alpha), DEC_LIMIT(1:n_alpha)
      double precision coord_illum(1:n_alpha, 1:n_delta)
      double precision day_night(1:n_alpha, 1:n_delta)
	integer j, k

c	Get the norm of the vector SAT to get its distance to Earth's center.
      call R_3D(SAT, R_SAT)

c     We look for the limit coordinates of the Earth seen by the satellite
c	alpha is the half-angle cone centred on the satellite containing the Earth
c	which is the horizon as seen by the satellite.

c	alpha is computed from the centre of the Earth !!
      alpha = acos(R_Earth/R_SAT)

c     We look for the coordinates of the Earth seen by the satellite

	do k = 1, n_alpha
	  alpha_grid = 2.d0*pi*dble(k-1)/dble(n_alpha-1)
	  do j = 1, n_delta
	    delta_grid = -pi/2.d0+dble(j-1)*pi/dble(n_delta-1)

	    dist = SphericalDistance(RA_SAT,DEC_SAT,alpha_grid,delta_grid)

	    if (dist <= alpha .and. day_night(k,j).ne.night) then
		coord_illum(k,j) = seen
	    else
c		Initialise the rest of the table.
		coord_illum(k,j) = night
	    endif
	  enddo
	enddo

	if (log_all_data) then
	  call write_seen(coord_illum)
	endif

      return
      end

c     ==========================================================================
      subroutine PHOTONS_AT_TELESCOPE(coord_illum, 				!Input
     &     RA_STAR, DEC_STAR, SAT,
     &     coord_photons, photons_hit, angle_photons, pst_photons,error)!Output

      implicit none
      include 'constants'
      include 'parameters'

      double precision coord_photons(1:n_alpha,1:n_delta)
      double precision coord_illum(1:n_alpha,1:n_delta)
      double precision angle_photons(1:n_alpha,1:n_delta)
      double precision pst_photons(1:n_alpha,1:n_delta)
      double precision prod, d, alpha_grid, delta_grid,d_cross
      double precision V(1:3), SAT(1:3), CELL(1:3)
      double precision angle_dat, pst_result
	double precision RA_STAR, DEC_STAR, R_STAR, STAR(1:3), vprod(1:3)
      double precision angle(n_ang), pst(n_ang)
      double precision min_angle_pst
	double precision a,b

      integer i, j, k
      logical photons_hit, error

      common /PST/ angle, pst, min_angle_pst

      error = .false.

      photons_hit = .false.

      R_STAR = 1.d0

	call CARTESIAN_COORDINATES(RA_STAR, DEC_STAR, R_STAR, !Input
     &     star) !Output

      do i = 1, n_alpha

         do j = 1, n_delta

            if(coord_illum(i,j).eq.illum) then 
               alpha_grid= 2.d0*pi*dble(i-1)/dble(n_alpha-1)
               delta_grid= -pi/2.d0 + dble(j-1)*pi/dble(n_delta-1)

               call CARTESIAN_COORDINATES(alpha_grid, delta_grid,
     &              R_Earth ,   !Input
     &              CELL)	  !Output

c		   Get the normalised vector pointing to the cell from the sat.
               call VEC_SUBSTRACTION(SAT, CELL, V)
               call R_3D(V, d)
		   V = V / d

               call scalar_product(V, STAR, prod)
		   call vect_product(STAR, V, vprod)
               call R_3D(vprod, d_cross)

               if(datan(d_cross/prod).gt.0) then 
                  coord_photons(i,j) = night !no photons hit the telescope
               else
                 photons_hit = .true.
                 coord_photons(i,j) = hit !photons hit the telescope

                 a = datan(d_cross/prod)*RAD
		     b = datan2(d_cross,prod)*RAD

                 a = pi/2.+datan(d_cross/prod)
		     b = pi/2.-(datan2(d_cross,prod)-pi/2.)

		     angle_photons(i,j) = b

		     if(include_pst) then
			 angle_dat = angle_photons(i,j)

c			 Check that the point given can actually be observed according
c			 to a simple rule: the angle between the line of sight and the
c			 Earth cannot be less than the stray light exclusion angle.
                   if(angle_dat*RAD.lt. SL_exclusion_angle .or. 
     & angle_dat*RAD .gt. 90d0) then
			   error = .true.
                   endif

			 if(angle_dat.lt.min_angle_pst) then
			   write(*,*)'Angle photons< min angle pst',  
     &                 angle_dat*RAD, min_angle_pst*RAD
                     error = .true.
                   endif

c			 Get the value of the PST for that particular LOS angle                   
			 call calcul_pst(angle_dat, pst_result)


                   if(pst_result.lt.pst_limit) then
                     write(*,*) 'PST < pst_limit=', pst_limit
                     write(*,*) 'angle_photons', angle_dat*RAD
                     write(*,*) 'pst_result',pst_result
                     error = .true.
                   endif

                   pst_photons(i,j) = pst_result 
                 endif
               endif

            endif

         enddo
      enddo

	if (log_all_data .and. error) then
	  call write_contamination(coord_photons)
	  call write_angles(angle_photons)
	  call write_star(RA_STAR, DEC_STAR)
	endif

	if (log_non_nominal .and. error) then
	  call write_error_targets(RA_STAR, DEC_STAR)
	endif 

      return
      end

c     ==========================================================================
      subroutine TOTAL_PHOTONS_AT_TELESCOPE(coord_photons,angle_photons,
     &     pst_photons, SAT,RA_SAT,DEC_SAT, day_night, !Input
     &     straylight_flux)                            !Output

      implicit none

      include 'constants'
      include 'parameters'

      double precision coord_photons(1:n_alpha,1:n_delta),
     &      angle_photons(1:n_alpha,1:n_delta), 
     &      pst_photons(1:n_alpha,1:n_delta), 
     &      day_night(1:n_alpha,1:n_delta)
      double precision dist(1:n_alpha,1:n_delta),
     &     EARTH(1:3), V(1:3)
      double precision d, SAT(1:3), R_SAT
      double precision energy, straylight_flux
      double precision solar_ref, solar_constant_cm
      double precision alpha_grid, delta_grid
      double precision Delta_alpha, Delta_delta
      double precision surf
      double precision RA_SAT, DEC_SAT
	double precision SphericalDistance
	double precision angle
	double precision proj_surf
      integer i,j

	double precision prod

	Delta_alpha = 2.d0*pi/dble(n_alpha-1)
	Delta_delta = pi/dble(n_delta-1)
	surf = R_Earth**2.d0 * Delta_alpha * Delta_delta

	energy = 0.d0
	straylight_flux = 0.d0

      call R_3D(SAT, R_SAT)

      solar_constant_cm = solar_constant/M2 !maximum energy that hit Earth per cm2
      solar_ref = solar_constant_cm * albedo/pi !reflected energy per cm2
    
      do i = 1, n_alpha
         do j = 1, n_delta

            if(coord_photons(i,j).eq.hit) then

               alpha_grid = 2.d0*pi*dble(i-1)/dble(n_alpha-1)
               delta_grid = -pi/2.d0 + dble(j-1)*pi/dble(n_delta-1)

               call CARTESIAN_COORDINATES(alpha_grid, delta_grid,
     &              R_Earth ,   ! Input
     &              EARTH)      ! Output
               
               call VEC_SUBSTRACTION(SAT, EARTH, V)

               call R_3D(V, d)

               call scalar_product(V, EARTH, prod) 

c		   Get the surface area of the projected surface
		   angle = SphericalDistance(RA_SAT, DEC_SAT,
     &		 alpha_grid, delta_grid)

c	The first angle is angle between cell and the satellite
c	The second angle is the angle between the LOS and the cell
c	The 3rd is for the fact towards high delta we have smaller surface
		   proj_surf = surf*dcos(angle)*abs(prod/(R_Earth*d))
     &		 *dcos(delta_grid)

               if(include_pst) then
                 energy = energy + 
     &                 solar_ref*proj_surf*pst_photons(i,j)*
     &                 day_night(i,j)/d**2.d0
               else
                  energy = energy + 
     &                 solar_ref*proj_surf*angle_photons(i,j)*
     &                 day_night(i,j)/d**2.d0
               endif       

            endif
         enddo
      enddo

      straylight_flux = energy * pixel !photons se px-1

      return
      end

c     =========================================================================
      subroutine vect_product(X, Y, Z)

      implicit none
      double precision X(1:3), Y(1:3), Z(1:3)
      

      Z(1) = X(2)*Y(3) - X(3)*Y(2)
      Z(2) = X(3)*Y(1) - X(1)*Y(3)
      Z(3) = X(1)*Y(2) - X(2)*Y(1)

      return
      end
c     ==========================================================================
      subroutine scalar_product(X, Y, z)

      implicit none
      double precision X(1:3), Y(1:3), z
      

      z = X(1)*Y(1) + X(2)*Y(2) + X(3)*Y(3)

      return
      end

c     ==========================================================================
      function REV(x)

      implicit none
      double precision REV, x

      REV = x - idint(x/360.d0)*360.d0
      if(REV.lt.0.d0) REV = REV + 360.d0

      return
      end

c     ==========================================================================
      function REV_RAD(x)

      implicit none
      double precision REV_RAD, x, twopi
      include 'constants'

      twopi = 2.d0*pi
      REV_RAD = x - idint(x/twopi)*twopi
      if(REV_RAD.lt.0.d0) REV_RAD = REV_RAD + twopi

      return
      end

c     =========================================================================
      subroutine RIGHT_ASCENSION(x, y, RA)

c     Given the cartesian coordinates it calculates the right ascention

      implicit none
      double precision x, y, RA
      
      RA = datan2(y,x)

      return
      end

c     =========================================================================
      subroutine DECLINATION(x, y, z, DEC)

c     Given the cartesian coordinates it calculates the declination

      implicit none
      double precision x, y, z, DEC

      DEC = datan2(z, sqrt(x**2+y**2))

      return
      end
c     =========================================================================   
      subroutine R_3D(X, R)

c     Given the cartesian coordinates it calculates the distance

      implicit none
      double precision X(1:3), R

      R = sqrt(X(1)**2+X(2)**2+X(3)**2)

      return
      end

c     =========================================================================   
      subroutine CARTESIAN_COORDINATES(alpha, delta, r, !I
     &     V)                                           !O

      implicit none
      double precision alpha, delta, r
      double precision V(1:3)

      V(1) = r*dcos(alpha)*dcos(delta)
      V(2) = r*dsin(alpha)*dcos(delta)
      V(3) = r*dsin(delta)

      return
      end
c     ========================================================================= 
      subroutine VEC_SUBSTRACTION(V, W, Z)

      implicit none
      double precision V(1:3), W(1:3), Z(1:3)

      Z(1) = V(1) - W(1)
      Z(2) = V(2) - W(2)
      Z(3) = V(3) - W(3)

      return
      end
      

c     =========================================================================   
      subroutine distance(R1, R2, D)

      implicit none
      double precision R1, R2, D

      D = sqrt(R1**2.d0 - R2**2.d0)

      return
      end

c     =========================================================================   
      subroutine calcul_pst(angle_dat, pst_result)

      implicit none
      include 'parameters'
      double precision angle_dat, pst_result

      call interp_lineal(angle_dat, pst_result)

      return
      end
c     =========================================================================  

      subroutine interp_lineal(angle_dat, pst_result)
 
      implicit none
      include 'parameters'
      double precision angle_dat, pst_result
      double precision angle(n_ang), pst(n_ang)
      double precision min_angle_pst
      common /PST/ angle, pst, min_angle_pst

      integer i

      do i= 1, n_ang-1
         if(angle_dat.ge.angle(i).and.angle_dat.lt.angle(i+1)) then 
            pst_result = pst(i) + 
     &     (angle_dat-angle(i))*(pst(i+1)-pst(i))/(angle(i+1)-angle(i))
         endif
      enddo

      return
      end


