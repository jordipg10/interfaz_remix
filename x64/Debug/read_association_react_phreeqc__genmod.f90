        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:53:12 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_ASSOCIATION_REACT_PHREEQC__genmod
          INTERFACE 
            SUBROUTINE READ_ASSOCIATION_REACT_PHREEQC(THIS,STRING,      &
     &PRIM_FLAG,DEFINED_SPECIES)
              USE EQ_REACTION_M, ONLY :                                 &
     &          EQ_REACTION_C,                                          &
     &          AQ_SPECIES_C
              CLASS (EQ_REACTION_C) :: THIS
              CHARACTER(*), INTENT(IN) :: STRING
              LOGICAL(KIND=4), INTENT(OUT) :: PRIM_FLAG
              TYPE (AQ_SPECIES_C) ,OPTIONAL, INTENT(OUT) ::             &
     &DEFINED_SPECIES
            END SUBROUTINE READ_ASSOCIATION_REACT_PHREEQC
          END INTERFACE 
        END MODULE READ_ASSOCIATION_REACT_PHREEQC__genmod
