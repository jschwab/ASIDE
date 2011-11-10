MODULE wh91

  USE types

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE accel(m1,m2,j1,j2)

    type(particle) :: m1,m2,j1,j2
    real(rl), dimension(3) :: dir1, indir1, dir2, indir2
    real(rl), dimension(3) :: r12

    j1 % a = 0d0

    r12 = (m2 % x - m1 % x)
    dir2 = - m1 % m * r12 / DOT_PRODUCT(r12,r12)**1.5d0

    indir2 = j2 % mu  / DOT_PRODUCT(j2 % x,j2 % x)**1.5d0 * &
             ((j2%x-m2%x) - f(j2%x-m2%x,m2%x) * m2%x)

    j2%a = dir2 + indir2

  END SUBROUTINE accel

  FUNCTION f(A,B)
    real(rl) :: f, q
    real(rl), dimension(3) :: A,B
    q = DOT_PRODUCT(A,A+2d0*B) / DOT_PRODUCT(B,B)
    f = q * (3d0 + q * (3d0 +q)) / (1d0 + (1d0 + q)**1.5d0)
  END FUNCTION F

END MODULE wh91
