        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:50 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRUESDELL_JONES_SUB__genmod
          INTERFACE 
            SUBROUTINE TRUESDELL_JONES_SUB(THIS,IONIC_ACT,LOG_ACT_COEFF)
              USE AQ_SPECIES_M
              CLASS (AQ_SPECIES_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: IONIC_ACT
              REAL(KIND=8), INTENT(OUT) :: LOG_ACT_COEFF
            END SUBROUTINE TRUESDELL_JONES_SUB
          END INTERFACE 
        END MODULE TRUESDELL_JONES_SUB__genmod
