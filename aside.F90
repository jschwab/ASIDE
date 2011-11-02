PROGRAM aside

  USE const
  USE types
  USE drift
  USE particles
  USE kick

  IMPLICIT NONE

  TYPE(particle) :: m0, m1, m2, j1, j2
  INTEGER ::  j
  REAL(rl) :: dt, a,e,i,o,w,f

  ! intialize perturbing particle
  m1 % m  = 1d-3
  m1 % mu = 1d0
  
  a = 1d0; e = 0d0; i = 0d0
  o = 0d0; w = 0d0; f = pi

  CALL oe_to_xv(m1,a,e,i,o,w,f)
 ! write(*,*) m1 %x, m1 %v

  ! intialize test mass
  m2 % m  = 0d0
  m2 % mu = 1d0
  
  a = 1d-1; e = 0d0; i = 15d0 * d2r
  o = 0d0; w = 0d0; f = 0d0

  CALL oe_to_xv(m2,a,e,i,o,w,f)
!  write(*,*) m2 %x, m2 %v

  CALL tojacobi(m1,m2,j1,j2)

  dt = 1d-1
  do j = 1, 100000
     call xv_to_oe(m2, a,e,i,o,w,f)
     write(*,"(I06, 9ES12.3)") j, a, e, i, o, w, f, m2%x
     call drift_kep(dt, j1)
     call drift_kep(dt, j2)
     call tohelio(m1,m2,j1,j2)
     call kick_all(dt, m1,m2,j1,j2)

  end do

  
END PROGRAM aside
