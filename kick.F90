MODULE kick

  USE types

  IMPLICIT NONE

  CONTAINS

!=====================================================================    

    SUBROUTINE kick_all(dt,m1,m2,j1,j2)

      USE wh91

      real(rl) :: dt
      type(particle) :: m1,m2,j1,j2

      CALL a(m1,m2,j1,j2)

      j1 % v = (j1 % v) + j1%a * dt
      j2 % v = (j2 % v) + j2%a * dt

    END SUBROUTINE kick_all

!=====================================================================    
  
END MODULE kick
