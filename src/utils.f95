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

!=======================================================================

  FUNCTION Cjacobi(p)
    
    TYPE(allparticles), INTENT(IN) :: p
    REAL(rl) :: Cjacobi
    REAL(rl) :: R0,R1,V2,n
    REAL(rl), dimension(3) :: H

    ! calculate the jacobi constant
    R0 = SQRT(DOT_PRODUCT(p%m2%x,p%m2%x))
    R1 = SQRT(DOT_PRODUCT(p%m1%x-p%m2%x,p%m1%x-p%m2%x))
    V2 = DOT_PRODUCT(p%m2%v,p%m2%v)
    
    n = SQRT(p%j1 % mu)

    H = cross(p%j2%x,p%j2%v)

    Cjacobi = -0.5d0*(v2 - 2d0 * (1d0/R0 + p%m1%m/R1) - 2d0 * H(3)*n)

   END FUNCTION Cjacobi


END MODULE utils
