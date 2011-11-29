MODULE wh91

  USE types

  IMPLICIT NONE

  CONTAINS

    SUBROUTINE accel(p)

      TYPE(ALLPARTICLES), INTENT(INOUT) :: p

      REAL(rl), DIMENSION(3) :: dir1, indir1
      REAL(rl), DIMENSION(3) :: dir2, indir2

      REAL(rl), DIMENSION(3) :: r12, r01, r02
      REAL(rl) :: r12m3

    ! calculate interparticle vectors r_ij -> r_j - r_i

      r12 = (p%m2%x - p%m1%x)
      r02 = (p%m2%x - p%m0%x)
      r01 = (p%m1%x - p%m0%x)

    ! pre-calculate magnitude(s) cubed

      r12m3 = DOT_PRODUCT(r12,r12)**(-1.5d0)

    ! calculate direct terms

      dir1 =   p%m2%m * r12 * r12m3
      dir2 = - p%m1%m * r12 * r12m3 & 
             - p%m1%m/p%j1%sigma * dir1 

    ! calculate indirect kicks

      indir1 = p%j1%mu * ( &
                  ((p%j1%x - r01) - f(p%j1%x,r01) * r01) &
                  / DOT_PRODUCT(p%j1 % x,p%j1 % x)**1.5d0 &
                 + p%m2%m/p%j1%sigma * r02 / DOT_PRODUCT(r02,r02)**1.5d0)

      indir2 = p%j2%mu * ( & 
                 ((p%j2%x - r02) - f(p%j2%x,r02) * r02) & 
                 / DOT_PRODUCT(p%j2 % x,p%j2 % x)**1.5d0 )

    ! set kicks

      p%j0%a = 0
      p%j1%a = dir1 + indir1 
      p%j2%a = dir2 + indir2

    end subroutine accel

  ! this function lets us rewrite the indirect terms such that
  ! we are unlikely subtract two approximately equal quantities
  ! see for example Battin (1987) pg. 389

    FUNCTION f(xj,xm)
      REAL(rl) :: f, q
      REAL(rl), DIMENSION(3), INTENT(IN) :: xj,xm
      q = DOT_PRODUCT(xj-xm,xj+xm) / DOT_PRODUCT(xm,xm)
      f = q * (3d0 + q * (3d0 +q)) / (1d0 + (1d0 + q)**1.5d0)
    END FUNCTION f

END MODULE wh91
