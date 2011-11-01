MODULE types

  IMPLICIT NONE

  REAL*8 :: real8
  INTEGER, PARAMETER :: rl = KIND(real8)

  TYPE particle
     REAL(rl) :: m  ! particle mass
     REAL(rl) :: mu  ! grav constant
     REAL(rl), DIMENSION(3) :: x, v ! jacobi phase space
     REAL(rl), DIMENSION(3) :: a ! kick
  END TYPE particle

END MODULE types
