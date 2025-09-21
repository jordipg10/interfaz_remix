!> Immobile zone module
module imm_zone_m
    implicit none
    save
    type, public :: imm_zone_c !> immobile zone class
        real(kind=8) :: imm_por !> immobile porosity (phi)
        real(kind=8) :: exch_rate !> exchange rate (alpha)
        real(kind=8) :: res_time !> residence time (tau)
        real(kind=8) :: prob !> probability
        !real(kind=8) :: flux !> flux of immobile zone
        !real(kind=8) :: conc
    end type
    
!> PFLOTRAN:
    
    !type, public :: immobile_type
    !
    !>    PetscInt :: nimmobile
    !>    PetscBool :: print_all
    !
    !>    type(immobile_species_type), pointer :: list
    !>    type(immobile_decay_rxn_type), pointer :: decay_rxn_list
    !
    !>    !> immobile species
    !>    character(len=MAXWORDLENGTH), pointer :: names(:)
    !>    PetscBool, pointer :: print_me(:)
    !
    !>    !> decay rxn
    !>    PetscInt :: ndecay_rxn
    !>    PetscInt, pointer :: decayspecid(:)
    !>    PetscReal, pointer :: decay_rate_constant(:)
    !end type immobile_type
end module