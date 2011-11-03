MODULE utils

  USE types
  
  IMPLICIT NONE
     
  CONTAINS

  FUNCTION cross(a, b) 

    REAL(rl), DIMENSION(3), INTENT(IN) :: a, b
    REAL(rl), DIMENSION(3) :: cross

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)

  END FUNCTION cross

END MODULE utils
