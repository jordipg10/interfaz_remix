        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:52:52 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_B_CONC_MOB__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_B_CONC_MOB(THIS,THETA,DELTA_T,           &
     &CONC_MOB_OLD,CONC_IMM_OLD,B)
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE PDE_MODEL_M
              USE MRMT_M, ONLY :                                        &
     &          MRMT_C
              CLASS (MRMT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(IN) :: CONC_MOB_OLD(:)
              REAL(KIND=8), INTENT(IN) :: CONC_IMM_OLD(:)
              REAL(KIND=8), INTENT(OUT) :: B(:)
            END SUBROUTINE COMPUTE_B_CONC_MOB
          END INTERFACE 
        END MODULE COMPUTE_B_CONC_MOB__genmod
