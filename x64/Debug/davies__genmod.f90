        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:23:57 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DAVIES__genmod
          INTERFACE 
            SUBROUTINE DAVIES(THIS,IONIC_ACT,A,LOG_ACT_COEFF)
              USE AQ_SPECIES_M, ONLY :                                  &
     &          AQ_SPECIES_C
              CLASS (AQ_SPECIES_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: IONIC_ACT
              REAL(KIND=8), INTENT(IN) :: A
              REAL(KIND=8), INTENT(OUT) :: LOG_ACT_COEFF
            END SUBROUTINE DAVIES
          END INTERFACE 
        END MODULE DAVIES__genmod
