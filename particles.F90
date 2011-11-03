MODULE particles
  
  USE types
  USE const
  USE utils
 
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE tojacobi(m1,m2,j1,j2)
 
    TYPE(particle), intent(in) :: m1, m2
    TYPE(particle), intent(inout) :: j1, j2
    REAL(rl) :: s0, s1, s2

    s0 = 1d0
    s1 = m1 % m + s0
    s2 = m2 % m + s1

    j1 % m = s0 / s1 * m1 % m
    j2 % m = s1 / s2 * m2 % m

    j1 % x = m1 % x
    j1 % v = m1 % v

    j2 % x = m2 % x - 1d0/s1 * (m1 % m * m1 % x)
    j2 % v = m2 % v - 1d0/s1 * (m1 % m * m1 % v)
    
    j1 % mu = s1 / s0
    j2 % mu = s2 / s1

  END SUBROUTINE tojacobi

!=======================================================================

  SUBROUTINE tohelio(m1,m2,j1,j2)
 
    TYPE(particle), intent(inout) :: m1, m2
    TYPE(particle), intent(in) :: j1, j2

    REAL(rl) :: s0,s1,s2

    s0 = 1d0
    s1 = m1 % m + s0
    s2 = m2 % m + s1

    m1 % x = j1 % x
    m1 % v = j1 % v

    m2 % x = j2 % x + m1 % m / s1 * j1 % x
    m2 % v = j2 % v + m1 % m / s1 * j1 % v
    
  END SUBROUTINE tohelio

!=======================================================================

  SUBROUTINE xv_to_oe(p,a,e,i,o,w,f)

  ! convert x,v to the standard orbital elements
  ! this uses the algorithm on page 53 of Murray & Dermott

  ! a,e,i,o,w,f - orbital elements
  ! p - particle whose x,v will be set
 
    TYPE(particle), INTENT(IN) :: p
    REAL(rl),INTENT(OUT):: a,e,i,o,w,f

    REAL(rl) :: R2,R,V2,H2,H,Rdot,wmf
    REAL(rl), DIMENSION(3) :: Hv

    R2 = DOT_PRODUCT(p%x,p%x)
    R  = SQRT(R2)
    V2 = DOT_PRODUCT(p%v,p%v)
    Hv = cross(p%x, p%v)
    H2 = DOT_PRODUCT(Hv,Hv)
    H  = SQRT(H2)

    Rdot = DOT_PRODUCT(p%x,p%v) / R

    a = 1d0 / (2d0/R - V2 / p%mu)
    e = SQRT(ABS(1d0 - H2/p%mu/a)) ! this is a cheap trick
    i = ACOS(Hv(3)/H)

    o = ATAN2(-Hv(1),Hv(2))

    wmf = ATAN2(p%x(3)*COS(o), p%x(1)*SIN(i)+p%x(3)*SIN(o)*COS(i))
    
    f = ATAN2(R*Rdot/H, 1d0 - R/a/(1d0-e*e))
    
    w = wmf + f
    
  END SUBROUTINE xv_to_oe

!=======================================================================
  
  SUBROUTINE oe_to_xv(p, a,e,i,o,w,f)

  ! convert the standard orbital elements to x,v
  ! this uses various expressions from  Murray & Dermott
  
  ! a,e,i,o,w,f - orbital elements
  ! p - particle whose x,v will be set

    TYPE(particle), INTENT(INOUT) :: p 
    REAL(rl), INTENT(IN) :: a,e,i,o,w,f

    REAL(rl), DIMENSION(3) :: x,v

    REAL(rl) :: ci,si,co,so,cwf,swf, naome2
    REAL(rl) :: r,rdot, rfdot

    ci = COS(i)
    si = SIN(i)
    co  = COS(o)
    so  = SIN(o)
    cwf = COS(w+f)
    swf = SIN(w+f)

    naome2 = SQRT(p % mu / a / (1d0 - e*e))

    x(1) = co * cwf - so * swf * ci
    x(2) = so * cwf + co * swf * ci
    x(3) = swf * si

    v(1) = -co * swf - so * cwf * ci
    v(2) = -so * swf + co * cwf * ci
    v(3) = cwf * si

    r = a * (1d0 - e*e) / (1d0 + e * COS(f))
    rdot = naome2 * e * SIN(f) 
    rfdot = naome2 * (1d0 + e*COS(f))
    
    p % v = rdot * x + rfdot * v
    p % x = r * x

  END SUBROUTINE oe_to_xv

!=======================================================================

END MODULE particles
