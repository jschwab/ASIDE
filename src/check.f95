MODULE check
       
  USE types

  IMPLICIT NONE

  CONTAINS

!=======================================================================

  SUBROUTINE check_close(p, tooclose)

    ! p  - particles to drift
    ! tooclose - are we too close?

    TYPE(allparticles), INTENT(INOUT) :: p
    LOGICAL :: tooclose
    REAL(rl) :: rh2, dx2

    ! this is your definition of a "few"
    REAL(rl), PARAMETER :: f = 2.5

    tooclose = .FALSE.

    ! separation of the particles
    dx2 = DOT_PRODUCT(p%m1%x-p%m2%x,p%m1%x-p%m2%x)

    ! size of the hill sphere
    rh2 = ((p%m1%m+p%m2%m)/3.0_rl)**(2.0_rl/3.0_rl)
    
    ! We're too close!
    IF (dx2 .LT. (f*f*rh2)) tooclose = .TRUE.      

  END SUBROUTINE check_close

!=======================================================================

END MODULE check
