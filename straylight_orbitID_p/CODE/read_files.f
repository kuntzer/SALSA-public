      subroutine READ_PST   

      implicit none
      include 'constants'
      include 'parameters'

      double precision  angle(1:n_ang), pst(1:n_ang)
      double precision  min_angle_pst
      common /PST/ angle, pst, min_angle_pst
      integer i

      open(unit= 20, file = pst_file_path)
      do i = 1, n_ang
         read(20, *) angle(i), pst(i)
      enddo
      close(20)
      
      do i = 1, n_ang
         angle(i) = angle(i)/RAD
      enddo

      min_angle_pst = angle(1)

      return
      END
c-----------------------------------------------------------------------
      subroutine READ_DATE_IN_JD(JD)

      implicit none
      double precision JD

      open(unit = 20, file= '../INPUT/date_JD.dat')
      read(20,*) JD
      close(20)

      return
      END
c-----------------------------------------------------------------------
      subroutine READ_COORD_SATELLITE(COORD_SAT, RA_SAT, DEC_SAT, R_SAT) 

      implicit none
      include 'constants'
      double precision COORD_SAT(1:3), RA_SAT, DEC_SAT, R_SAT

      open(unit= 20, file= '../INPUT/coord_sat.dat')

      read(20,*) COORD_SAT(1)
      read(20,*) COORD_SAT(2)
      read(20,*) COORD_SAT(3)
      read(20,*) RA_SAT
      read(20,*) DEC_SAT
      read(20,*) R_SAT

      close(20)

      COORD_SAT(1) = COORD_SAT(1)*KM
      COORD_SAT(2) = COORD_SAT(2)*KM
      COORD_SAT(3) = COORD_SAT(3)*KM
      R_SAT =  R_SAT*KM

      return
      END
c-----------------------------------------------------------------------
      subroutine READ_NUMBER_TARGETS(N_TARGETS)
 
      implicit none
      integer N_TARGETS

      open(unit = 20, file= '../INPUT/number_of_targets.dat')
      read(20,*) N_TARGETS
      close(20)

      return
      END
c-----------------------------------------------------------------------
      subroutine READ_TARGET_LIST(N_TARGETS, ALFA, DELTA)

      implicit none
      include 'parameters'
      double precision ALFA(n_targets_max), DELTA(n_targets_max)
      integer i, N_TARGETS

	if (N_TARGETS.gt.n_targets_max) then
		print *, "ERROR: The RA,DEC table is not big enough"
	endif

      do i = 1, N_TARGETS
         ALFA(i) = -20.d0
         DELTA(i) = -20.d0
      enddo
 
      open(unit = 20, file= '../INPUT/coord_targets.dat')

      do i = 1, N_TARGETS
      read(20,*) ALFA(i), DELTA(i) 
      enddo

      close(20)

c 100  format (f5.3, 1x, f5.3)

      return
      END
c-----------------------------------------------------------------------
      subroutine READ_COORD_SUN(RA_SUN, DEC_SUN, COORD_SUN) 

      implicit none
      include 'constants'
      double precision COORD_SUN(1:3), RA_SUN, DEC_SUN, R_SUN

      open(unit= 21, file= '../INPUT/coord_sun.dat')

      read(21,*) COORD_SUN(1)
      read(21,*) COORD_SUN(2)
      read(21,*) COORD_SUN(3)
      read(21,*) RA_SUN
      read(21,*) DEC_SUN
      read(21,*) R_SUN 

      COORD_SUN(1) = COORD_SUN(1)*KM
      COORD_SUN(2) = COORD_SUN(2)*KM
      COORD_SUN(3) = COORD_SUN(3)*KM

      close(21)

      return
      END
c-----------------------------------------------------------------------
