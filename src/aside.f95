PROGRAM aside

  USE const
  USE types
  USE utils
  USE particles
  USE step
  USE corrector
  USE output


  IMPLICIT NONE

  TYPE(allparticles) :: p
  TYPE(particle) :: m1, m2, j1, j2
  INTEGER ::  j
  REAL(rl) :: dt, a,e,i,o,w,f

  INTEGER, PARAMETER :: outfile = 7
  
  REAL(rl) :: R0,R1,V2,aa,n,CJ
  REAL(rl), dimension(3) :: H, xcm, vcm

  ! intialize perturbing particle
  m1 % m  = 1d-3
  m1 % mu = 1d0
  
  a = 1d0; e = 0d0; i = 0d0
  o = 0d0; w = 0d0; f =pi

  CALL oe_to_xv(m1,a,e,i,o,w,f)

  ! intialize test mass
  m2 % m  = 0d0
  m2 % mu = 1d0
  
  a = 0.1d0; e = 0.0d0; i = 30d0 * d2r
  o = 0d0; w = 0d0; f = 0.0d0

  CALL oe_to_xv(m2,a,e,i,o,w,f)

  CALL tojacobi(m1,m2,j1,j2)

  p % m1 = m1
  p % j1 = j1
  p % m2 = m2
  p % j2 = j2

  ! main loop

  open(unit=outfile, file = "aside.out", action = "write")

  CALL write_state(outfile, 0d0, p)

  dt = 1d-2
  do j = 1, 10000

     CALL C(dt,p)
     CALL bigstep(100,dt,p)
     CALL Ci(dt,p)
     CALL write_state(outfile, dt*j*100, p)
     
  end do

  close(outfile)
  
END PROGRAM aside


