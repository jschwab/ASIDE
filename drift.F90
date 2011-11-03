MODULE drift
       
  USE types

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE drift_kep(dt, p)

    USE kepler

    ! dt - integrator timestep
    ! p  - particle to drift

    REAL(rl), INTENT(IN) :: dt
    TYPE(particle), INTENT(INOUT) :: p
    
    REAL(rl), DIMENSION(3) :: x, v
    REAL(rl) :: f, g, fdot, gdot

    CALL gaussfg(dt, p, f, g, fdot, gdot)

    x = (p % x) * f    + (p % v) * g
    v = (p % x) * fdot + (p % v) * gdot

    p % x = x
    p % v = v

  END SUBROUTINE drift_kep

END MODULE drift
