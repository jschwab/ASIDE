MODULE step

  USE types
  USE coord
  USE drift
  USE kick

  IMPLICIT NONE

  CONTAINS

    SUBROUTINE bigstep(n, dt, p)

      ! n - number of substeps to take
      ! dt - timestep
      ! p - particles

      INTEGER, INTENT(IN) :: n 
      REAL(rl), INTENT(IN) :: dt
      TYPE(allparticles), INTENT(INOUT) :: p
      
      INTEGER :: i

      ! first part of step: D^{1/2}
      CALL drift_kep(0.5d0*dt, p)
      CALL tohelio(p)       

      DO i = 1, n - 1
         
         ! middle part of step (KD)^{n-1}

         CALL kick_all(dt, p)
         CALL drift_kep(dt,p)
         CALL tohelio(p)

      ENDDO

      ! final part of step KD^{1/2}
      CALL kick_all(dt, p)
      CALL drift_kep(0.5d0*dt, p)
      CALL tohelio(p)

    END SUBROUTINE bigstep

END MODULE step
