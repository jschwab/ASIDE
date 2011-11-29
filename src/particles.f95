MODULE particles
  
  USE types
  USE const
  USE utils
 
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE tojacobi(p)
 
    TYPE(allparticles), intent(inout) :: p

  ! this is called eta in M&D
    p%j0%sigma = 1.0d0
    p%j1%sigma = 1.0d0 + p%m1%m
    p%j2%sigma = 1.0d0 + p%m1%m + p%m2%m

  ! blah blah
    p%j0%m = 1d0 + p%m1%m + p%m2%m
    p%j1%m = p%m1%m / (1d0 + p%m1%m)
    p%j2%m = (1d0 + p%m1%m) * p%m2%m / (1d0 + p%m1%m + p%m2%m)

    p%j0%mu = 0d0
    p%j1%mu = 1d0 + p%m1%m
    p%j2%mu = (1d0 + p%m1%m + p%m2%m) / (1d0 + p%m1%m)

    p%j0%x = (p%m0%x + p%m1%m * p%m1%x + p%m2%m*p%m2%x)/(1d0+p%m1%m+p%m2%m)
    p%j0%v = (p%m0%v + p%m1%m * p%m1%v + p%m2%m*p%m2%v)/(1d0+p%m1%m+p%m2%m)
          
    p%j1%x = (p%m1%x-p%m0%x) 
    p%j1%v = (p%m1%v-p%m0%v) 
          
    p%j2%x = p%m2%x - (p%m0%x + p%m1%m * p%m1%x)/(1d0+p%m1%m)
    p%j2%v = p%m2%v - (p%m0%v + p%m1%m * p%m1%v)/(1d0+p%m1%m)


  END SUBROUTINE tojacobi

!=======================================================================

  SUBROUTINE tohelio(p)
 
    TYPE(allparticles), INTENT(INOUT) :: p

    p%m0%x = p%j0%x - p%m1%m * p%j1%x/(1d0+p%m1%m)-p%m2%m*p%j2%x/(1d0+p%m1%m+p%m2%m)
    p%m0%v = p%j0%v - p%m1%m * p%j1%v/(1d0+p%m1%m)-p%m2%m*p%j2%v/(1d0+p%m1%m+p%m2%m)

    p%m1%x = p%j0%x + p%j1%x/(1d0+p%m1%m)-p%m2%m*p%j2%x/(1d0+p%m1%m+p%m2%m)
    p%m1%v = p%j0%v + p%j1%v/(1d0+p%m1%m)-p%m2%m*p%j2%v/(1d0+p%m1%m+p%m2%m)

    p%m2%x = p%j0%x + (1d0 + p%m1%m) / (1d0 + p%m1%m + p%m2%m) * p%j2%x
    p%m2%v = p%j0%v + (1d0 + p%m1%m) / (1d0 + p%m1%m + p%m2%m) * p%j2%v

  END SUBROUTINE tohelio

!=======================================================================

  SUBROUTINE xv_to_oe(p,a,e,i,o,w,f)

  ! convert x,v to the standard orbital elements
  ! this uses the algorithm on page 53 of Murray & Dermott

  ! a,e,i,o,w,f - orbital elements
  ! p - particle whose x,v will be set
 
    TYPE(mparticle), INTENT(IN) :: p
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

    a = 1d0 / (2d0/R - V2)
    e = SQRT(MAX(1d0 - H2/a, 0d0))
       
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

    TYPE(mparticle), INTENT(INOUT) :: p 
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

    naome2 = SQRT(1d0 / a / (1d0 - e*e))

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

  FUNCTION Cjacobi(p)
    
    TYPE(allparticles), INTENT(IN) :: p
    REAL(rl) :: Cjacobi
    REAL(rl) :: R0,R1,V2,n
    REAL(rl), dimension(3) :: H

    ! calculate the jacobi constant
    R0 = SQRT(DOT_PRODUCT(p%m2%x,p%m2%x))
    R1 = SQRT(DOT_PRODUCT(p%m1%x-p%m2%x,p%m1%x-p%m2%x))
    V2 = DOT_PRODUCT(p%j2%v,p%j2%v)
    
    n = SQRT(p%j1 % mu)

    H = cross(p%j2%x,p%j2%v)

    Cjacobi = -0.5d0*(v2 - 2d0 * (1d0/R0 + p%m1%m/R1) - 2d0 * H(3)*n)

   END FUNCTION Cjacobi

!=======================================================================

END MODULE particles
