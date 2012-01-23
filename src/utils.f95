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
    REAL(rl) :: Rc,Rt,Rn,V2,n, Vn2, an
    REAL(rl), dimension(3) :: R, V, H

    R = p%m1%x - p%j0%x
    V = p%m1%v - p%j0%v

    V2 = DOT_PRODUCT(V,V)
    H = cross(R,V)

    ! calculate the jacobi constant
    Rc = SQRT(DOT_PRODUCT(p%m1%x-p%m0%x,p%m1%x-p%m0%x))
    Rt = SQRT(DOT_PRODUCT(p%m1%x-p%m2%x,p%m1%x-p%m2%x))

    Rn = SQRT(DOT_PRODUCT(p%m2%x-p%m0%x,p%m2%x-p%m0%x))
    Vn2 = DOT_PRODUCT(p%m2%v-p%m0%v,p%m2%v-p%m0%v)
    an = 1.0_rl / (2.0_rl / Rn - Vn2 / (1.0_rl + p%m2%m))
    n = SQRT((1.0_rl + p%m2%m) / an**3)
    
    Cjacobi = 2.0_rl * (1.0_rl / Rc + p%m2%m/Rt) + &
              2.0_rl * n * (R(1)*V(2)-R(2)*V(1)) - V2
              
   END FUNCTION Cjacobi


END MODULE utils
