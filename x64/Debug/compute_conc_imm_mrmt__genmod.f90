        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:20 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_CONC_IMM_MRMT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_CONC_IMM_MRMT(THIS,THETA,CONC_IMM_OLD,   &
     &CONC_MOB_OLD,CONC_MOB_NEW,DELTA_T,CONC_IMM_NEW)
              USE BCS_M
              USE PDE_MODEL_M
              USE MRMT_M, ONLY :                                        &
     &          MRMT_C
              CLASS (MRMT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(IN) :: CONC_IMM_OLD(:)
              REAL(KIND=8), INTENT(IN) :: CONC_MOB_OLD(:)
              REAL(KIND=8), INTENT(IN) :: CONC_MOB_NEW(:)
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(OUT) :: CONC_IMM_NEW(:)
            END SUBROUTINE COMPUTE_CONC_IMM_MRMT
          END INTERFACE 
        END MODULE COMPUTE_CONC_IMM_MRMT__genmod
