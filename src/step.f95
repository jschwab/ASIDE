MODULE step

  USE types
  USE coord
  USE drift
  USE kick
  USE check

  IMPLICIT NONE

  CONTAINS

    SUBROUTINE bigstep(n, dt, p, badstep)

      ! n - number of substeps to take
      ! dt - timestep
      ! p - particles

      INTEGER, INTENT(IN) :: n 
      REAL(rl), INTENT(IN) :: dt
      TYPE(allparticles), INTENT(INOUT) :: p
      LOGICAL, INTENT(OUT) :: badstep
      
      INTEGER :: i
      LOGICAL :: tooclose

      badstep = .FALSE.

!==================== DRIFT - KICK - DRIFT ====================

      ! ! first part of step: D^{1/2}
      ! CALL drift_kep(0.5d0*dt, p)
      ! CALL tohelio(p)       

      ! DO i = 1, n - 1
         
      !    ! middle part of step (KD)^{n-1}

      !    CALL kick_all(dt, p)
      !    CALL drift_kep(dt,p)
      !    CALL tohelio(p)

      !    CALL check_close(p,tooclose)
      !    if (tooclose) badstep = .TRUE.

      !    if (badstep) return

      ! ENDDO

      ! ! final part of step KD^{1/2}
      ! CALL kick_all(dt, p)
      ! CALL drift_kep(0.5d0*dt, p)
      ! CALL tohelio(p)

      ! CALL check_close(p,tooclose)
      ! if (tooclose) badstep = .TRUE.

!============================= END ===========================

!==================== KICK - DRIFT - KICK ====================

      ! first part of step: K^{1/2}

      CALL tohelio(p)       
      CALL kick_all(0.5_rl * dt, p)

      DO i = 1, n - 1
         
         ! middle part of step (DK)^{n-1}

         CALL drift_kep(dt,p)
         CALL tohelio(p)
         CALL check_close(p,tooclose)

         if (tooclose) badstep = .TRUE.
         if (badstep) return

         CALL kick_all(dt, p)

      ENDDO

      ! final part of step KD^{1/2}
      CALL drift_kep(dt, p)
      CALL tohelio(p)
      CALL check_close(p,tooclose)

      if (tooclose) badstep = .TRUE.
      if (badstep) return

      CALL kick_all(0.5_rl * dt, p)

!============================= END ===========================

      return

    END SUBROUTINE bigstep

END MODULE step
