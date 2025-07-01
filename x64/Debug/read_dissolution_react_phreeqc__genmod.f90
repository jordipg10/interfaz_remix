        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:53:16 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_DISSOLUTION_REACT_PHREEQC__genmod
          INTERFACE 
            SUBROUTINE READ_DISSOLUTION_REACT_PHREEQC(THIS,STRING,PHASE)
              USE EQ_REACTION_M, ONLY :                                 &
     &          EQ_REACTION_C,                                          &
     &          PHASE_C
              CLASS (EQ_REACTION_C) :: THIS
              CHARACTER(*), INTENT(IN) :: STRING
              CLASS (PHASE_C) ,OPTIONAL, INTENT(INOUT) :: PHASE
            END SUBROUTINE READ_DISSOLUTION_REACT_PHREEQC
          END INTERFACE 
        END MODULE READ_DISSOLUTION_REACT_PHREEQC__genmod
