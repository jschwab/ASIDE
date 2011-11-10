MODULE corrector

  USE types
  USE drift
  USE kick

  IMPLICIT NONE

  REAL(rl), PARAMETER, PRIVATE :: alpha = 0.4183300132670377739890860d0
  REAL(rl), PARAMETER, PRIVATE :: beta = 0.04980119205559973499870072d0
  
  INTEGER, PARAMETER :: Cstage = 2
  REAL(rl), DIMENSION(2), PARAMETER, PRIVATE  :: a = (/1d0, -1d0/) * alpha
  REAL(rl), DIMENSION(2), PARAMETER, PRIVATE  :: b =  (/0.5d0, -0.5d0/) * beta

  CONTAINS

!=======================================================================

  SUBROUTINE C(dt,p)

    REAL(rl), INTENT(IN) :: dt
    TYPE(allparticles), INTENT(INOUT) :: p

    INTEGER :: i

    DO i = Cstage, 1, -1
       CALL z(dt*a(i),dt*b(i),p)
    END DO

  END SUBROUTINE C

!=======================================================================

  SUBROUTINE Ci(dt,p)

    REAL(rl), INTENT(IN) :: dt
    TYPE(allparticles), INTENT(INOUT) :: p

    INTEGER :: i

    DO i = 1, Cstage
       CALL z(dt*a(i),-dt*b(i),p)
    END DO

  END SUBROUTINE Ci

!=======================================================================

  SUBROUTINE z(adt,bdt,p)

    REAL(rl), INTENT(IN) :: adt,bdt
    TYPE(allparticles), INTENT(INOUT) :: p

    CALL chi(-adt,-bdt, p)
    CALL chi( adt, bdt, p)

  END SUBROUTINE z
   
!=======================================================================
 
  SUBROUTINE chi(adt,bdt,p)

    REAL(rl), INTENT(IN) :: adt,bdt
    TYPE(allparticles), INTENT(INOUT) :: p

    call drift_kep(adt,p)
    call kick_all(bdt,p)
    call drift_kep(adt,p)

  END SUBROUTINE chi

!=======================================================================

END MODULE corrector
