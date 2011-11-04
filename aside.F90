PROGRAM aside

  USE const
  USE types
  USE drift
  USE particles
  USE kick
  USE utils

  IMPLICIT NONE

  TYPE(particle) :: m0, m1, m2, j1, j2
  INTEGER ::  j
  REAL(rl) :: dt, a,e,i,o,w,f
  
  REAL(rl) :: R0,R1,V2,aa,n,CJ
  REAL(rl), dimension(3) :: H, xcm, vcm

  ! intialize perturbing particle
  m1 % m  = 1d-3
  m1 % mu = 1d0
  
  a = 1d0; e = 0d0; i = 0d0
  o = 0d0; w = 0d0; f =pi

  CALL oe_to_xv(m1,a,e,i,o,w,f)
 ! write(*,*) m1 %x, m1 %v

  ! intialize test mass
  m2 % m  = 0d0
  m2 % mu = 1d0
  
  a = 0.1d0; e = 0.0d0; i = 15d0 * d2r
  o = 0d0; w = 0d0; f = 0.0d0

  CALL oe_to_xv(m2,a,e,i,o,w,f)
!  write(*,*) m2 %x, m2 %v


  CALL tojacobi(m1,m2,j1,j2)

  dt = 1d-3
  do j = 1, 1000 * 100000

     if (MOD(j,1000) == 0) then 

     ! calculate the jacobi constant
     R0 = SQRT(DOT_PRODUCT(m2%x,m2%x))
     R1 = SQRT(DOT_PRODUCT(m1%x-m2%x,m1%x-m2%x))
     V2 = DOT_PRODUCT(j2%v,j2%v)

     n = SQRT(j1 % mu)

     H = cross(j2%x,j2%v)

     CJ = v2 - 2d0 * (1d0 / R0 + m1%m / R1) - 2d0 * H(3)*n

     call xv_to_oe(m2, a,e,i,o,w,f)
     write(*,"(11ES16.6)") j*dt, a, e, i, o, w, f, m2%x, CJ

     end if

     call drift_kep(0.5*dt, j1)
     call drift_kep(0.5*dt, j2)

     call tohelio(m1,m2,j1,j2)
     call kick_all(dt, m1,m2,j1,j2)

     call drift_kep(0.5*dt, j1)
     call drift_kep(0.5*dt, j2)



  end do

  
END PROGRAM aside
