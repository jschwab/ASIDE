MODULE kepler

  USE types
  USE const
  
  IMPLICIT NONE

  ! this module provides tools to solve kepler's equation
  ! to a specified tolerance using a maximum # of iterations

  ! TODO: make these user-editable
  REAL(rl), PARAMETER :: tol = 1e-12
  INTEGER, PARAMETER :: imax = 10

  CONTAINS

!=====================================================================

  SUBROUTINE danby4(x, deltaM, ecE0, esE0)

  ! solve difference form of Kepler's equation using an algorithm
  ! from Danby, J.A.M., "Fundamentals of Celestial Mechanics". 

  ! x - guess for deltaE -> final deltaE
  ! deltaM - delta mean anomaly (n * dt)
  ! ecE0 - e * cos(E0)
  ! esE0 - e * sin(E0)

    REAL(rl), INTENT(INOUT) :: x
    REAL(rl), INTENT(IN) :: deltaM, ecE0, esE0

    INTEGER :: i 
    REAL(rl) :: sx, cx
    REAL(rl) :: f, fp, fpp, fppp
    REAL(rl) :: dx

    DO i = 1, imax

       sx = SIN(x)
       cx = COS(x)

       f = x - ecE0*sx + esE0*(1d0 - cx) - deltaM
       fp = 1d0 - ecE0*cx + esE0*sx
       fpp = ecE0*sx + esE0*cx
       fppp = ecE0*cx - esE0*sx

       dx = -f / fp
       dx = -f /(fp + dx*fpp/2d0)
       dx = -f /(fp + dx*(fpp/2d0 + dx*fppp/6d0))
       
       x = x + dx
       IF (ABS(dx/x) <= tol) EXIT

    END DO

    RETURN

  END SUBROUTINE danby4

!=====================================================================

  SUBROUTINE gaussfg(dt, p, f, g, fdot, gdot)

  ! evaluate gauss f & g functions
 
  ! dt - timestep
  ! p - particle
  ! f,g - gauss f & g functions
  ! fdot, gdot - time derivatives of gauss f & g functions

    REAL(rl), INTENT(IN) :: dt
    TYPE(jparticle), INTENT(IN) :: p
    REAL(rl), INTENT(OUT) :: f, g, fdot, gdot

    REAL(rl) :: r0, v02, u, a, n, r, e, sigma
    REAL(rl) :: deltaE, deltaM, ecE0, esE0, cE, sE, y


    ! calculate the orbital elements
    r0 = SQRT(DOT_PRODUCT(p % x, p % x))
    v02 = DOT_PRODUCT(p % v, p % v)
    u = DOT_PRODUCT(p % x, p % v)
    a = 1d0 / (2d0/r0 - v02/p % mu)
    n = SQRT(p % mu / a**3)

    deltaM = n * dt
    ecE0 = 1d0 - r0/a
    esE0 = u/(n*a*a)

    e = SQRT(ecE0*ecE0 + esE0*esE0)

    ! danby suggests the following guess for delta E
    IF (e .lt. 0.1) THEN
       deltaE = deltaM
    ELSE ! high eccentricity
       y = MOD(deltaM,twopi) - esE0
       sigma = SIGN(0.85d0, esE0*COS(y) + ecE0*SIN(y))
       deltaE = y + sigma * e
    END IF

    ! calculate delta E
    CALL danby4(deltaE, deltaM, ecE0, esE0)
  
    cE = COS(deltaE)
    sE = SIN(deltaE)
    
    ! TODO: make the computation of f & g more efficient
    ! r 
    r = a * (1d0 - ecE0 * cE + esE0 * sE)

    ! gauss f & g
    f = (a/r0) * (cE - 1d0) + 1d0
    g = dt + (1d0/n) * (sE - deltaE)

    ! gauss fdot & gdot
    fdot = -(a/r)*(a/r0)* n * sE
    gdot = (a/r)*(cE - 1d0) + 1
    
    RETURN
    
  END SUBROUTINE gaussfg

!=====================================================================

END MODULE kepler
