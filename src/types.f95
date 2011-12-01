MODULE types

  IMPLICIT NONE

  ! rl = 8 -> double precision
  ! rl = 16 -> quad precision (requires libquadmath, for example)

  INTEGER, PARAMETER :: rl = 8

  ! "real" particles, in heliocentric coordinates
  TYPE mparticle
     REAL(rl) :: m  ! particle mass
     REAL(rl), DIMENSION(3) :: x, v ! phase space
  END TYPE mparticle

  ! "effective" particles, in Jacobi coordinates
  TYPE jparticle
     REAL(rl) :: m  ! particle mass
     REAL(rl) :: sigma  ! particle mass
     REAL(rl) :: mu  ! particle mass
     REAL(rl), DIMENSION(3) :: x, v ! phase space
     REAL(rl), DIMENSION(3) :: a ! kick
  END TYPE jparticle

  ! container for all the particles
  TYPE allparticles
     TYPE(mparticle) :: m0 ! central mass
     TYPE(mparticle) :: m1
     TYPE(mparticle) :: m2
     TYPE(jparticle) :: j0 ! center of mass 
     TYPE(jparticle) :: j1 
     TYPE(jparticle) :: j2 
  END TYPE allparticles

END MODULE types
