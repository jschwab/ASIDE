PROGRAM aside

  USE const
  USE types
  USE utils
  USE step
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

  CALL write_state(outfile, 0.0_rl, p)

  ndump = iNT(tdump / dt)
  nmax  = INT(tmax / dt)

  t = 0

  DO n = 1, nmax, ndump

     CALL bigstep(ndump,dt,p)

     t = t + dt*dble(ndump)

     CALL write_state(outfile, t, p)
     
  END DO

  CLOSE(outfile)
  
END PROGRAM aside


