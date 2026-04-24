!> \file lin_kin_reaction_m.f90
!> \brief Linear kinetic reaction module.
!> \details
!> Defines the `lin_kin_reaction_c` class for first-order (linear) kinetic reactions
!> of the form \f$ A \to B \f$ with rate law:
!> \f[
!>   r_k(c) = \lambda \, c
!> \f]
!> where \f$\lambda\f$ [1/T] is the first-order rate constant and \f$c\f$ [M/L^3]
!> is the concentration of the reacting species.
!>
!> The Jacobian is simply \f$ \partial r_k / \partial c = \lambda \f$.
!>
!> \see kin_reaction_m, redox_kin_reaction_m
!> \author Jordi
!> \date Unknown
!> \ingroup kinetics
module lin_kin_reaction_m
    use kin_reaction_m, only : kin_reaction_c
    use aq_phase_m, only : aq_phase_c
    use species_m, only : species_c
    implicit none
    save
    private
    !> \brief Linear kinetic reaction subclass.
    !> \details First-order decay/growth reaction with rate constant \f$\lambda\f$.
    type, public, extends(kin_reaction_c) :: lin_kin_reaction_c
        real(kind=8) :: lambda  !< [1/T] First-order rate constant
    contains
        procedure :: set_lambda                !< Set rate constant \f$\lambda\f$
        procedure :: set_index_aq_phase_lin     !< Set aqueous phase index for reacting species
        procedure :: compute_rk_lin            !< Compute \f$ r_k = \lambda c \f$
        procedure :: compute_drk_dc_lin        !< Compute \f$ \partial r_k / \partial c = \lambda \f$
        procedure :: write_params=>write_lambda !< Write rate constant to file
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
                
        subroutine set_index_aq_phase_lin(this,aq_phase,species)
            implicit none
            class(lin_kin_reaction_c) :: this
            class(aq_phase_c), intent(in) :: aq_phase
            class(species_c), intent(in) :: species(:)
            
            integer(kind=4) :: aq_species_ind
            logical :: flag
            
            allocate(THIS%indices_aq_phase(1))
            call aq_phase%is_species_in_aq_phase(species(this%species_ind(1)),flag,aq_species_ind)
            if (flag.eqv..true.) then
                this%indices_aq_phase(1)=aq_species_ind
            else
                error stop "Linear species is not in aqueous phase"
            end if
        end subroutine
        
        subroutine write_lambda(this,unit)
            implicit none
            class(lin_kin_reaction_c), intent(in) :: this
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