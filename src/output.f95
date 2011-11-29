MODULE OUTPUT

  USE types
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
      TYPE(allparticles), INTENT(IN) :: p

      REAL(rl), DIMENSION(3) :: x, v
      REAL(rl), DIMENSION(6) :: oe

      x = p%m2%x - p%m0%x
      v = p%m2%v - p%m0%v

      CALL xv_to_oe(x,v,oe)
!      CJ = Cjacobi(p)

      WRITE(outfile,"(7ES20.10)") t, p%m2%x-p%m0%x,p%m2%v-p%m0%v
!      WRITE(outfile,"(7ES20.10)") t, oe
      
    END SUBROUTINE write_state

!=======================================================================

END MODULE OUTPUT
