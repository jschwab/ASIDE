MODULE types

  IMPLICIT NONE

  INTEGER, PARAMETER :: rl = 8

  TYPE mparticle
     REAL(rl) :: m  ! particle mass
     REAL(rl), DIMENSION(3) :: x, v ! jacobi phase space
  END TYPE mparticle

  TYPE jparticle
     REAL(rl) :: m  ! particle mass
     REAL(rl) :: sigma  ! particle mass
     REAL(rl) :: mu  ! particle mass
     REAL(rl), DIMENSION(3) :: x, v ! jacobi phase space
     REAL(rl), DIMENSION(3) :: a ! kick
  END TYPE jparticle


  TYPE allparticles
     TYPE(mparticle) :: m0 ! central mass
     TYPE(mparticle) :: m1 ! perturbing mass
     TYPE(mparticle) :: m2 ! test mass
     TYPE(jparticle) :: j0 
     TYPE(jparticle) :: j1 
     TYPE(jparticle) :: j2 

  END TYPE allparticles

END MODULE types
