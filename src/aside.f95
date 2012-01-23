PROGRAM aside

  USE const
  USE types
  USE utils
  USE step
  USE input
  USE output

  IMPLICIT NONE

  TYPE(allparticles) :: p
  LOGICAL :: badstep

  INTEGER ::  n
  REAL(rl) :: t, dt, tmax, tdump

  INTEGER :: nmax, ndump

  INTEGER, PARAMETER ::  infile = 7
  INTEGER, PARAMETER :: outfile = 8

  ! read the initial conditions
  OPEN(unit=infile, file = "aside.in", action = "read")
  CALL read_input(infile, p, dt, tmax, tdump)
  CLOSE(infile)

  ! open the output file and dump the intial state
  OPEN(unit=outfile, file = "aside.out", action = "write")
  CALL write_state(outfile, 0.0_rl, p)

  ! convert the time intervals into step intervals
  ndump = INT(tdump / dt)
  nmax  = INT(tmax / dt)

  ! we start at t = 0 and then off we go
  t = 0
  DO n = 1, nmax, ndump

     ! we advance  ndump steps at once
     ! this is only approximately tdump
     t = t + dt*REAL(ndump,rl)
     CALL bigstep(ndump,dt,p, badstep)
     WRITE(6,"(A, ES8.2, A, I3,A)") "t = ", t, " (", INT(t/tmax * 1e2_rl),"%)"

     ! if something bad happened, it's time to stop
     IF (badstep) THEN 
        WRITE(6,*) "Terminating..."
        EXIT
     ENDIF

     ! we want output, so here it is
     CALL write_state(outfile, t, p)
     
  END DO

  ! and we're all finished now
  CLOSE(outfile)
  
END PROGRAM aside


