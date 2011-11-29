MODULE coord
  
  USE types
  USE const
  USE utils
 
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE tojacobi(p)
 
    TYPE(allparticles), intent(inout) :: p

  ! this is called eta in M&D
    p%j0%sigma = 1.0d0
    p%j1%sigma = 1.0d0 + p%m1%m
    p%j2%sigma = 1.0d0 + p%m1%m + p%m2%m

  ! blah blah
    p%j0%m = 1d0 + p%m1%m + p%m2%m
    p%j1%m = p%m1%m / (1d0 + p%m1%m)
    p%j2%m = (1d0 + p%m1%m) * p%m2%m / (1d0 + p%m1%m + p%m2%m)

    p%j0%mu = 0d0
    p%j1%mu = 1d0 + p%m1%m
    p%j2%mu = (1d0 + p%m1%m + p%m2%m) / (1d0 + p%m1%m)

    p%j0%x = (p%m0%x + p%m1%m * p%m1%x + p%m2%m*p%m2%x)/(1d0+p%m1%m+p%m2%m)
    p%j0%v = (p%m0%v + p%m1%m * p%m1%v + p%m2%m*p%m2%v)/(1d0+p%m1%m+p%m2%m)
          
    p%j1%x = (p%m1%x-p%m0%x) 
    p%j1%v = (p%m1%v-p%m0%v) 
          
    p%j2%x = p%m2%x - (p%m0%x + p%m1%m * p%m1%x)/(1d0+p%m1%m)
    p%j2%v = p%m2%v - (p%m0%v + p%m1%m * p%m1%v)/(1d0+p%m1%m)


  END SUBROUTINE tojacobi

!=======================================================================

  SUBROUTINE tohelio(p)
 
    TYPE(allparticles), INTENT(INOUT) :: p

    p%m0%x = p%j0%x - p%m1%m * p%j1%x/(1d0+p%m1%m)-p%m2%m*p%j2%x/(1d0+p%m1%m+p%m2%m)
    p%m0%v = p%j0%v - p%m1%m * p%j1%v/(1d0+p%m1%m)-p%m2%m*p%j2%v/(1d0+p%m1%m+p%m2%m)

    p%m1%x = p%j0%x + p%j1%x/(1d0+p%m1%m)-p%m2%m*p%j2%x/(1d0+p%m1%m+p%m2%m)
    p%m1%v = p%j0%v + p%j1%v/(1d0+p%m1%m)-p%m2%m*p%j2%v/(1d0+p%m1%m+p%m2%m)

    p%m2%x = p%j0%x + (1d0 + p%m1%m) / (1d0 + p%m1%m + p%m2%m) * p%j2%x
    p%m2%v = p%j0%v + (1d0 + p%m1%m) / (1d0 + p%m1%m + p%m2%m) * p%j2%v

  END SUBROUTINE tohelio

!=======================================================================

END MODULE coord
