!> Kinetic reaction module
module kin_reaction_m
    use reaction_m, only: reaction_c, write_reaction_sup
    use kin_params_m, only: kin_params_c
    implicit none
    save
    type, public, abstract, extends(reaction_c) :: kin_reaction_c !> Kinetic reaction abstract subclass
        integer(kind=4) :: num_aq_rk=0 !> number of aqueous species relevant for kinetic reaction rate
        integer(kind=4), allocatable :: indices_aq_phase(:) !> indices of species in aqueous phase object that control kinetic reaction rates
    contains
        procedure, public :: set_num_aq_rk
        procedure, public :: allocate_indices_aq_phase_kin_react
        procedure(write_params), public, deferred :: write_params
    end type kin_reaction_c
    
    type, public :: kin_reaction_poly_c !< clase ad hoc para crear vector punteros clase reaccion cinetica (cosas de Fortran)
        class(kin_reaction_c), pointer :: kin_reaction
    contains
        procedure, public :: set_kin_reaction
    end type kin_reaction_poly_c
    
    abstract interface
        subroutine write_params(this,unit)
            import kin_reaction_c
            implicit none
            class(kin_reaction_c) :: this
            integer(kind=4), intent(in) :: unit !> file unit
        end subroutine
        
    end interface
    
    interface 
        
        
        
    end interface
    
    contains
        subroutine set_kin_reaction(this,kin_reaction)
            class(kin_reaction_poly_c) :: this
            class(kin_reaction_c), intent(in), target :: kin_reaction
            this%kin_reaction=>kin_reaction
        end subroutine

        subroutine set_num_aq_rk(this,num_aq_rk)
            class(kin_reaction_c) :: this
            integer(kind=4), intent(in) :: num_aq_rk
            if (num_aq_rk<0) then
                error stop "Number of aqueous species relevant for kinetic reaction rate must be positive"
            else
                this%num_aq_rk=num_aq_rk
            end if
        end subroutine

        subroutine allocate_indices_aq_phase_kin_react(this,num_aq_rk)
            class(kin_reaction_c) :: this
            integer(kind=4), intent(in), optional :: num_aq_rk
            if (present(num_aq_rk)) then
                call this%set_num_aq_rk(num_aq_rk)
            end if
            allocate(this%indices_aq_phase(this%num_aq_rk))
        end subroutine
            
end module kin_reaction_m