MODULE OUTPUT

  USE types
  USE coord
  USE orbel

  IMPLICIT NONE

  CONTAINS

!=======================================================================

    SUBROUTINE write_state(outfile, t,p)
      
      ! outfile - 
      ! t - time
      ! p - particles

      INTEGER, INTENT(IN) :: outfile
      REAL(rl), INTENT(IN) :: t
      TYPE(allparticles), INTENT(INOUT) :: p

      REAL(rl), DIMENSION(3) :: x, v
      REAL(rl), DIMENSION(6) :: oe

      REAL(rl) :: CJ

      CALL tohelio(p)

      x = p%m1%x - p%m0%x
      v = p%m1%v - p%m0%v

      CALL xv_to_oe(x,v,oe)
      CJ = Cjacobi(p)

!      WRITE(outfile,"(7ES20.10)") t, x, v
      WRITE(outfile,"(7ES20.10)") t, oe(1:5), CJ
      
    END SUBROUTINE write_state

!=======================================================================

END MODULE OUTPUT
