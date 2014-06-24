c=======================================================================

      subroutine write_illuminated_earth(day_night)

      implicit none
      include 'parameters'
      include 'constants'
      double precision alpha_grid, delta_grid
      double precision day_night(1:n_alpha, 1:n_delta)
      integer i,j 



      open(unit = 21, file = '../OUTPUT/illum_earth.out')

      do i = 1, n_alpha
         alpha_grid = 2.d0*pi*dble(i-1)/dble(n_alpha-1)
         
         do j = 1, n_delta
            delta_grid = -pi/2.d0 + dble(j-1)*pi/dble(n_delta-1)

            write(21,*) alpha_grid*RAD, delta_grid*RAD, 
     &           day_night(i,j)*albedo 
            write(21,*)
         enddo
      enddo

      close(21)

      return
      end

c=======================================================================

      subroutine write_seen(day_night)

      implicit none
      include 'parameters'
      include 'constants'
      double precision alpha_grid, delta_grid
      double precision day_night(1:n_alpha, 1:n_delta)
      integer i,j 



      open(unit = 21, file = '../OUTPUT/seen.out')

      do i = 1, n_alpha
         alpha_grid = 2.d0*pi*dble(i-1)/dble(n_alpha-1)
         
         do j = 1, n_delta
            delta_grid = -pi/2.d0 + dble(j-1)*pi/dble(n_delta-1)

            write(21,*) alpha_grid*RAD, delta_grid*RAD, 
     &           day_night(i,j)*albedo 
            write(21,*)
         enddo
      enddo

      close(21)

      return
      end

c=======================================================================

      subroutine write_angles(day_night)

      implicit none
      include 'parameters'
      include 'constants'
      double precision alpha_grid, delta_grid
      double precision day_night(1:n_alpha, 1:n_delta)
      integer i,j 



      open(unit = 21, file = '../OUTPUT/angles.out')

      do i = 1, n_alpha
         alpha_grid = 2.d0*pi*dble(i-1)/dble(n_alpha-1)
         
         do j = 1, n_delta
            delta_grid = -pi/2.d0 + dble(j-1)*pi/dble(n_delta-1)

            write(21,*) alpha_grid*RAD, delta_grid*RAD, 
     &           day_night(i,j)*albedo 
            write(21,*)
         enddo
      enddo

      close(21)

      return
      end

c=======================================================================

      subroutine write_contamination(day_night)

      implicit none
      include 'parameters'
      include 'constants'
      double precision alpha_grid, delta_grid
      double precision day_night(1:n_alpha, 1:n_delta)
      integer i,j 



      open(unit = 21, file = '../OUTPUT/contamination.out')

      do i = 1, n_alpha
         alpha_grid = 2.d0*pi*dble(i-1)/dble(n_alpha-1)
         
         do j = 1, n_delta
            delta_grid = -pi/2.d0 + dble(j-1)*pi/dble(n_delta-1)

            write(21,*) alpha_grid*RAD, delta_grid*RAD, 
     &           day_night(i,j)*albedo 
            write(21,*)
         enddo
      enddo

      close(21)

      return
      end

c=======================================================================

      subroutine write_star(ra,dec)

      implicit none
      include 'parameters'
      include 'constants'
      double precision ra,dec

      open(unit = 21, file = '../OUTPUT/star.out')

      write(21,*) ra, dec

      close(21)

      return
      end

c=======================================================================

      subroutine write_error_targets(ra,dec)

      implicit none
      include 'parameters'
      include 'constants'
      double precision ra,dec, JD

      call READ_DATE_IN_JD(JD)

      open(unit = 21, file = '../OUTPUT/errors.out', ACCESS='APPEND')

      write(21,*) JD, ra, dec

      close(21)

      return
      end
