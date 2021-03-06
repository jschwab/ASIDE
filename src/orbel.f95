MODULE orbel
  
  USE types
  USE const
  USE utils
 
  IMPLICIT NONE

  REAL(rl), PARAMETER :: TINY = 1e-12_rl

  ! the orbital elements are manipulated in a 6-element array
  ! (a,e,i,o,w,f)
  ! a - semimajor axis
  ! e - eccentricity
  ! i - inclination
  ! o - longitude of ascending node
  ! w - argument of pericenter
  ! f - true anomaly

  CONTAINS

!=======================================================================

  SUBROUTINE xv_to_oe(mu,x,v,oe)

  ! convert x,v to the standard orbital elements
  ! this uses the algorithm on page 53 of Murray & Dermott
 
    REAL(rl), INTENT(IN) :: mu
    REAL(rl), DIMENSION(3), INTENT(IN) :: x, v
    REAL(rl), DIMENSION(6), INTENT(OUT) :: oe

    REAL(rl) :: a,e,i,o,w,f
    REAL(rl) :: R2,R,V2,H2,H,Rdot,wpf
    REAL(rl), DIMENSION(3) :: Hv

    REAL(rl) :: fac, cape, sw, cw

    R2 = DOT_PRODUCT(x,x)
    R  = SQRT(R2)
    V2 = DOT_PRODUCT(v,v)
    Hv = cross(x, v)
    H2 = DOT_PRODUCT(Hv,Hv)
    H  = SQRT(H2)

    Rdot = DOT_PRODUCT(x,v) / R

    ! semi-major axis
    a = 1.0_rl / (2.0_rl / R - V2/mu)

    ! eccentricity
    fac = 1d0 - H2/a/mu
    IF (fac > TINY) THEN
       e = SQRT(fac)
    ELSE
       e = 0.0_rl
    ENDIF

    ! inclination
    i = ACOS(Hv(3)/H)

    ! longitude of ascending node
    fac = SQRT(Hv(1)*Hv(1)+Hv(2)*Hv(2))/H
    IF (fac > TINY) THEN
       o = ATAN2(Hv(1),-Hv(2))
    ELSE
       o = 0d0
    ENDIF

    ! put o in [0, 2 pi)
    IF (o .LT. 0.0_rl) o = o + twopi

    ! argument of pericenter
    wpf = ATAN2(x(3)*COS(o), x(1)*SIN(i)+x(3)*SIN(o)*COS(i))    

    ! true anomaly
    f = ATAN2(a * (1.0_rl - e*e) * Rdot / H, &
              a * (1.0_rl - e*e)/R - 1.0_rl)

    IF (o == 0d0) then 
       w = 0
    ELSE
       w = wpf - f
    END IF

    IF (w .LT. 0) w = w + twopi

    oe(1) = a
    oe(2) = e
    oe(3) = i
    oe(4) = o
    oe(5) = w
    oe(6) = f

    
  END SUBROUTINE xv_to_oe

!=======================================================================
  
  SUBROUTINE oe_to_xv(mu,x,v,oe)

  ! convert the standard orbital elements to x,v
  ! this uses various expressions from  Murray & Dermott
  
    REAL(rl), INTENT(IN) :: mu
    REAL(rl), DIMENSION(3), INTENT(OUT) :: x, v
    REAL(rl), DIMENSION(6), INTENT(IN) :: oe

    REAL(rl) :: a,e,i,o,w,f
    REAL(rl) :: ci,si,co,so,cwf,swf, naome2
    REAL(rl) :: r,rdot, rfdot

    a = oe(1) 
    e = oe(2) 
    i = oe(3) 
    o = oe(4) 
    w = oe(5) 
    f = oe(6) 

    ci = COS(i)
    si = SIN(i)
    co  = COS(o)
    so  = SIN(o)
    cwf = COS(w+f)
    swf = SIN(w+f)

    naome2 = SQRT(mu / a / (1d0 - e*e))

    x(1) = co * cwf - so * swf * ci
    x(2) = so * cwf + co * swf * ci
    x(3) = swf * si

    v(1) = -co * swf - so * cwf * ci
    v(2) = -so * swf + co * cwf * ci
    v(3) = cwf * si

    r = a * (1d0 - e*e) / (1d0 + e * COS(f))
    rdot = naome2 * e * SIN(f) 
    rfdot = naome2 * (1d0 + e*COS(f))
    
    v = rdot * x + rfdot * v
    x = r * x

  END SUBROUTINE oe_to_xv

!=======================================================================

 END MODULE orbel
