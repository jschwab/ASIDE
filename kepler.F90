module kepler

  ! this module provides tools to solve kepler's equation
  ! to a specified tolerance using a maximum # of iterations

  ! TODO: make these user-editable
  real(rl), parameter :: tol = 1e-12
  integer, parameter :: imax = 10

  subroutine danby4(x, deltaM, ecE0, esE0)

  ! solve difference form of Kepler's equation using the algorithm in
  ! Danby
  ! which has local quartic convergence

  ! x - guess for deltaE -> final deltaE
  ! deltaM - delta mean anomaly (n * dt)
  ! ecE0 - e * cos(E0)
  ! esE0 - e * sin(E0)

    real(rl), intent(inout) :: x
    real(rl), intent(in) :: deltaM,ecE0, esE0
    real(rl) :: sx, cx
    real(rl) :: f, fp, fpp, fppp
    real(rl) :: dx
    
    do i = 1, imax

       sx = sin(x)
       cx = cos(x)

       f = x - ecE0*sx + esE0*(1d0 - cx) - deltaM
       fp = 1d0 - ecE0*c + esE0*sx
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


end module kepler
