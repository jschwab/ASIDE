module kepler

  use types
  use const
  
  implicit none

  ! this module provides tools to solve kepler's equation
  ! to a specified tolerance using a maximum # of iterations

  ! TODO: make these user-editable
  real(rl), parameter :: tol = 1e-12
  integer, parameter :: imax = 10

  contains

  subroutine danby4(x, deltaM, ecE0, esE0)

  ! solve difference form of Kepler's equation using the algorithm in
  ! Danby
  ! which has local quartic convergence

  ! x - guess for deltaE -> final deltaE
  ! deltaM - delta mean anomaly (n * dt)
  ! ecE0 - e * cos(E0)
  ! esE0 - e * sin(E0)


    real(rl), intent(inout) :: x
    real(rl), intent(in) :: deltaM, ecE0, esE0

    integer :: i 
    real(rl) :: sx, cx
    real(rl) :: f, fp, fpp, fppp
    real(rl) :: dx

    do i = 1, imax

       sx = sin(x)
       cx = cos(x)

       f = x - ecE0*sx + esE0*(1d0 - cx) - deltaM
       fp = 1d0 - ecE0*cx + esE0*sx
       fpp = ecE0*sx + esE0*cx
       fppp = ecE0*cx - esE0*sx

       dx = -f / fp
       dx = -f /(fp + dx*fpp/2d0)
       dx = -f /(fp + dx*(fpp/2d0 + dx*fppp/6d0))
       
       x = x + dx
       if (abs(dx) <= tol) exit

    end do

    return

  end subroutine danby4


  subroutine gaussfg(dt, obj, f, g, fdot, gdot)

    use const
    use types

    real(rl), intent(in) :: dt
    type(particle), intent(in) :: obj
    real(rl), intent(out) :: f, g, fdot, gdot

    real(rl) :: r0, v02, u, a, n, r, e, sigma
    real(rl) :: deltaE, deltaM, ecE0, esE0, cE, sE, y


    ! calculate the orbital elements
    r0 = sqrt(dot_product(obj % x, obj % x))
    v02 = dot_product(obj % v, obj % v)
    u = dot_product(obj % x, obj % v)
    a = 1d0 / (2d0/r0 - v02/obj % mu)
    n = sqrt(obj % mu / a**3)

    deltaM = n * dt
    ecE0 = 1d0 - r0/a
    esE0 = u/(n*a*a)

    e = sqrt(ecE0*ecE0 + esE0*esE0)

    ! danby suggestes the following guess for delta E
    if (e .lt. 0.1) then
       deltaE = deltaM
    else ! high eccentricity
       y = mod(deltaM,twopi) - esE0
       sigma = sign(0.85d0, esE0*cos(y) + ecE0*sin(y))
       deltaE = y + sigma * e
    end if

    ! calculate delta E
    call danby4(deltaE, deltaM, ecE0, esE0)
  
    cE = cos(deltaE)
    sE = sin(deltaE)
    
    ! r 
    r = a * (1d0 - ecE0 * cE + esE0 * sE)

    ! gauss f & g
    f = (a/r0) * (cE - 1d0) + 1d0
    g = dt + (1d0/n) * (sE - deltaE)

    ! gauss fdot & gdot
    fdot = -(a/r)*(a/r0)* n * sE
    gdot = (a/r)*(cE - 1d0) + 1
    
    return
    
  end subroutine gaussfg


end module kepler
