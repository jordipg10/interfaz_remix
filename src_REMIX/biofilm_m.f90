!> Biofilm module (to be developed)
module biofilm_m
    use imm_zone_m
    use species_m
    implicit none
    save
    type, public, extends(imm_zone_c) :: biofilm_c !> biofilm class
        !> This class contains the information of the biofilm
        !> It extends the imm_zone class
        !> The immobile porosity is the volumetric fraction of the biofilm (CHAPUZA)
        real(kind=8) :: f_w !> volumetric fraction of water in the biofilm
        real(kind=8) :: yield !> yield factor (rate of biomass produced to substrate consumed) (OJO: YA LO TIENES EN LA CLASE REACCION REDOX)
        integer(kind=4) :: num_mcrg !> number of microorganisms
        type(microorganism_c), allocatable :: microorganisms(:) !> microorganisms
        real(kind=8), allocatable :: vol_fracs(:) !> microorganism volumetric fractions (QUIZAS MEJOR EN SOLID CHEMISTRY)
    contains
        procedure, public :: compute_biomass_decay
    end type
    
!> PFLOTRAN:
    
    !type, public :: immobile_species_type
    !>    PetscInt :: id
    !>    character(len=MAXWORDLENGTH) :: name
    !>    PetscReal :: molar_weight
    !>    PetscBool :: print_me
    !>    type(immobile_species_type), pointer :: next
    !end type immobile_species_type
    
    !type, public :: microbial_rxn_type
    !>    PetscInt :: id
    !>    PetscInt :: itype
    !>    character(len=MAXSTRINGLENGTH) :: reaction
    !>    PetscReal :: rate_constant
    !>    PetscReal :: activation_energy
    !>    PetscBool :: print_me
    !>    type(reaction_equation_type), pointer :: reaction_equation
    !>    type(monod_type), pointer :: monod
    !>    type(inhibition_type), pointer :: inhibition
    !>    type(microbial_biomass_type), pointer :: biomass
    !>    type(microbial_rxn_type), pointer :: next
    !>  end type microbial_rxn_type

    !type, public :: monod_type
    !>    PetscInt :: id
    !>    character(len=MAXWORDLENGTH) :: species_name
    !>    PetscReal :: half_saturation_constant
    !>    PetscReal :: threshold_concentration
    !>    type(monod_type), pointer :: next
    !>  end type monod_type
    
    !type, public :: microbial_biomass_type
    !>    PetscInt :: id
    !>    character(len=MAXWORDLENGTH) :: species_name
    !>    PetscReal :: yield
    !end type microbial_biomass_type
    
    contains
    
    function compute_biomass_decay(this,conc) result(rk_d)
        class(biofilm_c) :: this
        real(kind=8), intent(in) :: conc(:) !> concentration of microbes
        real(kind=8), allocatable :: rk_d(:) !> reaction rate of the biomass
        
        integer(kind=4) :: i
        
        allocate(rk_d(this%num_mcrg))
        do i=1,this%num_mcrg
            if (this%microorganisms(i)%decay_flag) then
                rk_d(i)=this%microorganisms(i)%decay_rate*conc(i)
            else
                rk_d(i)=0d0
            end if
        end do
    end function
end module