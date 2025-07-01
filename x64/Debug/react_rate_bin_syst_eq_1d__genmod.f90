        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:52:38 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE REACT_RATE_BIN_SYST_EQ_1D__genmod
          INTERFACE 
            SUBROUTINE REACT_RATE_BIN_SYST_EQ_1D(THIS,U,DU_DX,D,PHI,R_EQ&
     &)
              USE EQ_REACTION_M, ONLY :                                 &
     &          EQ_REACTION_C
              CLASS (EQ_REACTION_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: U
              REAL(KIND=8), INTENT(IN) :: DU_DX
              REAL(KIND=8), INTENT(IN) :: D
              REAL(KIND=8), INTENT(IN) :: PHI
              REAL(KIND=8), INTENT(OUT) :: R_EQ
            END SUBROUTINE REACT_RATE_BIN_SYST_EQ_1D
          END INTERFACE 
        END MODULE REACT_RATE_BIN_SYST_EQ_1D__genmod
