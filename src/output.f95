MODULE OUTPUT

  USE types
  USE particles

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

      REAL(rl) :: a,e,i,o,w,f,CJ

      CALL xv_to_oe(p%m2, a,e,i,o,w,f)
      CJ = Cjacobi(p)
      WRITE(outfile,"(11ES16.6)") t, a, e, i, o, w, f, p%m2%x, CJ

    END SUBROUTINE write_state

!=======================================================================

END MODULE OUTPUT
