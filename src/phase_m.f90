!> Phase superclass
!! This class contains the properties of a phase
module phase_m
    implicit none
    save
    type, public, abstract :: phase_c !> phase abstract superclass 
        integer(kind=4) :: num_species=0 !> number of species
        integer(kind=4) :: num_var_act_species=0 !> number of variable activity species
        integer(kind=4) :: num_cst_act_species=0 !> number of constant activity species
        character(len=256)  :: name !> phase name
    contains
    !> Set
        procedure, public :: set_num_species_phase
        procedure, public :: set_num_var_act_species_phase
        procedure, public :: set_num_cst_act_species_phase
        procedure, public :: set_phase_name
    !> Compute
        procedure, public :: compute_num_species_phase
    end type
    
    contains
        subroutine set_phase_name(this,name)
            implicit none
            class(phase_c) :: this
            character(len=*), intent(in) :: name
            this%name=name
        end subroutine 
        
        subroutine set_num_species_phase(this,num_species)
            implicit none
            class(phase_c) :: this
            integer(kind=4), intent(in) :: num_species
            if (num_species<0) error stop "Number of species cannot be negative"
            this%num_species=num_species
        end subroutine 
        
        subroutine compute_num_species_phase(this)
            implicit none
            class(phase_c) :: this
            this%num_species=this%num_var_act_species+this%num_cst_act_species
        end subroutine 
        
        subroutine set_num_var_act_species_phase(this,num_var_act_species)
            implicit none
            class(phase_c) :: this
            integer(kind=4), intent(in) :: num_var_act_species
            if (num_var_act_species<0) then
                error stop "Number of variable activity species cannot be negative"
            else if (num_var_act_species>this%num_species) then
                error stop "Number of variable activity species cannot be greater than number of species"
            else
                this%num_var_act_species=num_var_act_species
            end if
        end subroutine 
        
        subroutine set_num_cst_act_species_phase(this,num_cst_act_species)
            implicit none
            class(phase_c) :: this
            integer(kind=4), intent(in) :: num_cst_act_species
            if (num_cst_act_species<0) then
                error stop "Number of species cannot be negative"
            else if (num_cst_act_species>this%num_species) then
                error stop "Number of constant activity species cannot be greater than number of species"
            else
                this%num_cst_act_species=num_cst_act_species
            end if
        end subroutine 
end module