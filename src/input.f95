MODULE INPUT

  USE const
  USE types
  USE coord
  USE orbel

  IMPLICIT NONE

  CONTAINS

!=======================================================================

    SUBROUTINE read_input(infile, p,dt,tmax,tdump)

      INTEGER, INTENT(IN) :: infile
      REAL(rl), INTENT(OUT) :: dt,tmax,tdump
      TYPE(allparticles), INTENT(INOUT) :: p
      TYPE(mparticle) :: mc, mp, mt

      REAL(rl), DIMENSION(6) :: oe
      REAL(rl) :: m,a,e,i,o,w,f

      NAMELIST /pertmass/ m, a, e, i, o, w, f
      NAMELIST /testmass/ a, e, i, o, w, f
      NAMELIST /timestep/ dt, tmax, tdump

      ! timestep & output controls
      READ(infile,nml=timestep)
      
      ! units of mass are central mass
      ! the central mass is at the origin
      ! and it is instantaneously at rest
      mc % m = 1.0_rl
      mc % x = 0.0_rl
      mc % v = 0.0_rl

      ! get the information about the peturbing mass
      READ(infile,nml=pertmass)
      mp % m = m

      ! we read orbital elements; angles are in degrees
      oe(1) = a;       oe(2) = e;       oe(3) = i * d2r
      oe(4) = o * d2r; oe(5) = w * d2r; oe(6) = f * d2r

      CALL oe_to_xv(mp%x,mp%v,oe)

      ! get the information about the test mass
      READ(infile,nml=testmass)
      mt % m = 0_rl

      ! we read orbital elements; angles are in degrees
      oe(1) = a;       oe(2) = e;       oe(3) = i * d2r
      oe(4) = o * d2r; oe(5) = w * d2r; oe(6) = f * d2r

      CALL oe_to_xv(mt%x, mt%v, oe)

      mt % x = (/0.1_rl, 0.0_rl, 0.0_rl/)
      mt % v = (/0.0_rl, 3.0_rl, 0.0_rl/)

      mp % x = (/-1.0_rl, 0.0_rl, 0.0_rl/)
      mp % v = (/ 0.0_rl,-1.0_rl, 0.0_rl/)

      ! this chooses the ordering in jacobi coordinates
      p % m0 = mc
      p % m1 = mt
      p % m2 = mp

      ! convert to jacobi coordinatse
      CALL tojacobi(p)

    END SUBROUTINE read_input

!=======================================================================

END MODULE INPUT
