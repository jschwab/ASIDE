MODULE kick

  USE types

  IMPLICIT NONE

  CONTAINS

!=====================================================================    

    SUBROUTINE kick_all(dt,p)

      USE wh91

      REAL(rl),INTENT(IN) :: dt
      TYPE(allparticles),INTENT(INOUT) :: p

      CALL accel(p)
      
      p%j1 % v = (p%j1 % v) + p%j1%a * dt
      p%j2 % v = (p%j2 % v) + p%j2%a * dt

    END SUBROUTINE kick_all

!=====================================================================    
  
END MODULE kick
