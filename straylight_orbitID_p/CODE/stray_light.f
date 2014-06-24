      program STRAY_LIGHT

      implicit none
      include 'constants'
      include 'parameters'
      double precision COORD_SAT(1:3), RA_SAT, DEC_SAT, R_SAT
      double precision RA_SUN, DEC_SUN
      double precision RA_STAR, DEC_STAR
      double precision ALFA(1:n_targets_max), DELTA(1:n_targets_max)
	double precision copy_DELTA(1:n_targets_max)
      double precision coord(1:n_alpha, 1:n_delta),
     &     coord_illum(1:n_alpha, 1:n_delta), 
     &     coord_photons(1:n_alpha,1:n_delta),
     &     angle_photons(1:n_alpha,1:n_delta),
     &     pst_photons(1:n_alpha,1:n_delta),
     &     proj_surf(1:n_alpha,1:n_delta),
     &     day_night(1:n_alpha, 1:n_delta),
     &	   COORD_SUN(1:3)
      double precision straylight_flux

      integer N_TARGETS, i, count_error

      logical photons_hit, error

      photons_hit = .false.

c--------------------------------------------------------------------------------
c     READ INPUT FILES

      if(include_pst) call READ_PST

      call READ_COORD_SATELLITE(COORD_SAT,RA_SAT,DEC_SAT,R_SAT)

      call READ_NUMBER_TARGETS(N_TARGETS)

      call READ_TARGET_LIST(N_TARGETS, ALFA, DELTA)

      call READ_COORD_SUN(RA_SUN, DEC_SUN, COORD_SUN)
c--------------------------------------------------------------------------------
c     Coordinates of the illuminated part of the Earth by the SUN

      call ILLUMINATED_EARTH(RA_SUN, DEC_SUN, !I
     &     day_night)           !O

c--------------------------------------------------------------------------------
c     Calculate the coordinates of the Earth as seen from the satellite
      call EARTH_SEEN_FROM_SAT(COORD_SAT, RA_SAT, DEC_SAT, day_night, !Input
     &     coord_illum)               !Output

c--------------------------------------------------------------------------------
c     Calculate the coordinates of the illuminated Earth that the satellite sees

	count_error = 0

	open(unit=211, file='../OUTPUT/straylight.out')
      
	do i = 1, N_TARGETS

         DEC_STAR = DELTA(i)
         RA_STAR = ALFA(i)

c--------------------------------------------------------------------------------
c     Given the coordinates of the target star, we caculate the energy 
c     at the telescope for the given direction         

         call PHOTONS_AT_TELESCOPE(coord_illum, RA_STAR,	! Input
     &        DEC_STAR, COORD_SAT,
     &        coord_photons, photons_hit,				! Output
     & 	      angle_photons, pst_photons, error)

	if (.not. error) then 
         if(photons_hit) then
            
c     for the regions where the photons that hit the detector come we:
c     1) take into account the orientation of the satellite or include the PST

            call TOTAL_PHOTONS_AT_TELESCOPE(coord_photons,
     &           angle_photons,  pst_photons, COORD_SAT,
     &		RA_SAT,DEC_SAT, day_night,  !Input
     &           straylight_flux)          !Output
         else
            straylight_flux = 0.d0
         endif

c	2) Finally, we write the result for this particular direction into the output file.
	 write(211,*) RA_STAR,DEC_STAR,straylight_flux
	else 
		count_error = count_error + 1
	endif

      enddo

      close(211)! close flux file

	if (verbose .and. count_error .gt. 0) write(*,"(I2,a17)") 
     & count_error,' points discarded'

      END

