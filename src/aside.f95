PROGRAM aside

  USE const
  USE types
  USE utils
  USE particles
  USE step
  USE corrector
  USE input
  USE output


  IMPLICIT NONE

  TYPE(allparticles) :: p

  INTEGER ::  n
  REAL(rl) :: t, dt, tmax, tdump

  INTEGER :: nmax, ndump

  INTEGER, PARAMETER ::  infile = 7
  INTEGER, PARAMETER :: outfile = 8

  
  OPEN(unit=infile, file = "aside.in", action = "read")
  CALL read_input(infile, p, dt, tmax, tdump)
  CLOSE(infile)

  OPEN(unit=outfile, file = "aside.out", action = "write")

  CALL write_state(outfile, 0d0, p)

  ndump = int(tdump / dt)
  nmax  = int(ndump / dt)

  t = 0

  do n = 1, nmax

!     CALL C(dt,p)
     CALL bigstep(ndump,dt,p)
!     CALL Ci(dt,p)

     t = t + dt*dble(ndump)

     CALL write_state(outfile, t, p)
     
  end do

  close(outfile)
  
END PROGRAM aside


