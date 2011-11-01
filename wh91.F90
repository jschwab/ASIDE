MODULE wh91

  USE types

  CONTAINS

  SUBROUTINE a(m1,m2,j1,j2)

    type(particle) :: m1,m2,j1,j2
    real(rl), dimension(3) :: dir1, indir1, dir2, indir2
    real(rl), dimension(3) :: r12

    r12 = (m1 % x - m2 % x)
    dir2 = - m1 % m * r12 / DOT_PRODUCT(r12,r12)**1.5d0

    indir2 = - j2 % mu  / DOT_PRODUCT(j2 % x,j2 % x)**1.5 * &
            ((j2%x-m2%x) - f(j2%x-m2%x,m2%x) * m2%x)


    j1 % a = 0d0
    j2 % a = dir2 + indir2 
    
  END SUBROUTINE a

  FUNCTION f(r1,r2)
    real(rl) :: f, q
    real(rl), dimension(3) :: r1, r2
    q = DOT_PRODUCT(r1, r1 - 2d0 * r2) / DOT_PRODUCT(r2,r2)
    f = q * (3d0 + q * (3d0 +q)) / (1d0 + (1d0 + q)**1.5d0)
  END FUNCTION F

END MODULE wh91
