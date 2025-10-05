!> Linear kinetics module
!!> $A \to B$
!!> \f$r_k(c)=\lambda*c_\f$
module lin_kin_reaction_m
    use kin_reaction_m, only : kin_reaction_c
    use aq_phase_m, only : aq_phase_c
    implicit none
    save
    type, public, extends(kin_reaction_c) :: lin_kin_reaction_c !> linear kinetic reaction subclass
        real(kind=8) :: lambda !> reaction rate
    contains
    !> Set
        procedure, public :: set_lambda
        procedure, public :: set_index_aq_phase_lin
    !> Compute
        procedure, public :: compute_rk_lin
        procedure, public :: compute_drk_dc_lin
    !> Write
        procedure, public :: write_params=>write_lambda
    end type
    
    interface
                
    end interface
    
    contains
        subroutine set_lambda(this,lambda)
            implicit none
            class(lin_kin_reaction_c) :: this
            real(kind=8), intent(in) :: lambda
            this%lambda=lambda
        end subroutine
                
        subroutine set_index_aq_phase_lin(this,aq_phase)
            implicit none
            class(lin_kin_reaction_c) :: this
            class(aq_phase_c), intent(in) :: aq_phase
            
            integer(kind=4) :: aq_species_ind
            logical :: flag
            
            allocate(THIS%indices_aq_phase(1))
            call aq_phase%is_species_in_aq_phase(this%species(1),flag,aq_species_ind)
            if (flag.eqv..true.) then
                this%indices_aq_phase(1)=aq_species_ind
            else
                error stop "Linear species is not in aqueous phase"
            end if
        end subroutine
        
        subroutine write_lambda(this,unit)
            implicit none
            class(lin_kin_reaction_c) :: this
            integer(kind=4), intent(in) :: unit !> file unit
            write(unit,*) this%lambda
        end subroutine
        
        subroutine compute_drk_dc_lin(this,drk_dc)
            implicit none
            class(lin_kin_reaction_c), intent(in) :: this
            real(kind=8), intent(out) :: drk_dc(:)
            drk_dc=0d0
            drk_dc(this%indices_aq_phase)=this%lambda
        end subroutine
        
        subroutine compute_rk_lin(this,conc,rk)
            implicit none
            class(lin_kin_reaction_c), intent(in) :: this
            real(kind=8), intent(in) :: conc
            real(kind=8), intent(out) :: rk
            rk=this%lambda*conc
        end subroutine


        
end module