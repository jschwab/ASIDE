MODULE types

  IMPLICIT NONE

  INTEGER, PARAMETER :: rl = 8

  TYPE particle
     REAL(rl) :: m  ! particle mass
     REAL(rl) :: mu  ! grav constant
     REAL(rl), DIMENSION(3) :: x, v ! jacobi phase space
     REAL(rl), DIMENSION(3) :: a ! kick
  END TYPE particle

  TYPE allparticles
     TYPE(particle) :: m1, j1 ! perturbing mass
     TYPE(particle) :: m2, j2 ! test mass
  END TYPE allparticles

END MODULE types
