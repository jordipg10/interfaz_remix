        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:46 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_EQ_REACTION__genmod
          INTERFACE 
            SUBROUTINE READ_EQ_REACTION(THIS,SPECIES,FILENAME)
              USE PARAMS_SPEC_VOL_M
              USE PARAMS_ACT_COEFF_M
              USE SPECIES_M
              USE EQ_REACTION_M, ONLY :                                 &
     &          EQ_REACTION_C,                                          &
     &          SPECIES_C
              CLASS (EQ_REACTION_C) :: THIS
              CLASS (SPECIES_C), INTENT(IN) :: SPECIES
              CHARACTER(*), INTENT(IN) :: FILENAME
            END SUBROUTINE READ_EQ_REACTION
          END INTERFACE 
        END MODULE READ_EQ_REACTION__genmod
