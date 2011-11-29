MODULE drift
       
  USE types

  IMPLICIT NONE

  CONTAINS

!=======================================================================

  SUBROUTINE drift_kep(dt, p)

    ! dt - integrator timestep
    ! p  - particles to drift

    REAL(rl), INTENT(IN) :: dt
    TYPE(allparticles), INTENT(INOUT) :: p

    p % j0 % x = p % j0 % x + p % j0 % v * dt
    CALL drift_one(dt,p%j1)
    CALL drift_one(dt,p%j2)

  END SUBROUTINE drift_kep

!=======================================================================

  SUBROUTINE drift_one(dt, pn)

    USE kepler

    ! dt - integrator timestep
    ! p  - particles to drift

    REAL(rl), INTENT(IN) :: dt
    TYPE(jparticle), INTENT(INOUT) :: pn
    
    REAL(rl), DIMENSION(3) :: x, v
    REAL(rl) :: f, g, fdot, gdot

    CALL gaussfg(dt, pn, f, g, fdot, gdot)

    x = (pn % x) * f    + (pn % v) * g
    v = (pn % x) * fdot + (pn % v) * gdot

    pn % x = x
    pn % v = v

  END SUBROUTINE drift_one

!=======================================================================

END MODULE drift
