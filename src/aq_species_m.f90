!> Aqueous species subclass:
!!>   contains the properties of aqueous species
!!>   computes logarithm activity coefficient with different models
module aq_species_m
    use species_m
    implicit none
    save
    type, public, extends(species_c) :: aq_species_c !> aqueous species subclass
        real(kind=8) :: ionic_radius !> expressed in angstroms
        real(kind=8) :: alk_contrib !> contribution to alkalinity
    contains
    !> Models to compute activity coefficient
        procedure, public :: Debye_Huckel_restr
        procedure, public :: Debye_Huckel_ampl
        procedure, public :: Davies
        procedure, public :: Truesdell_Jones
    end type
    
    contains
        
        
        subroutine Debye_Huckel_restr(this,ionic_act,A,B,log_act_coeff)
            implicit none
            class(aq_species_c) :: this
            real(kind=8), intent(in) :: ionic_act !> I
            REAL(kind=8), intent(in) :: A,B
            real(kind=8), intent(out) :: log_act_coeff !> log_10(gamma)
            
            log_act_coeff=-A*sqrt(ionic_act)*this%valence**2
        end subroutine
        
        subroutine Debye_Huckel_ampl(this,ionic_act,A,B,log_act_coeff)
            implicit none
            class(aq_species_c) :: this
            real(kind=8), intent(in) :: ionic_act !> I
            REAL(kind=8), intent(in) :: A,B
            real(kind=8), intent(out) :: log_act_coeff !> log_10(gamma)
            
            log_act_coeff=-(A*sqrt(ionic_act)*this%valence**2)/(1d0+this%params_act_coeff%ion_size_param*B*sqrt(ionic_act))
        end subroutine
        
        subroutine Davies(this,ionic_act,A,B,log_act_coeff)
            implicit none
            class(aq_species_c) :: this
            real(kind=8), intent(in) :: ionic_act !> I
            REAL(kind=8), intent(in) :: A,B
            real(kind=8), intent(out) :: log_act_coeff !> log_10(gamma)
            
            log_act_coeff=-A*(this%valence**2)*((sqrt(ionic_act)/(1d0+sqrt(ionic_act)))-3d-1*ionic_act)
        end subroutine
        
        subroutine Truesdell_Jones(this,ionic_act,A,B,log_act_coeff)
            implicit none
            class(aq_species_c) :: this
            real(kind=8), intent(in) :: ionic_act !> I
            REAL(kind=8), intent(in) :: A,B
            real(kind=8), intent(out) :: log_act_coeff !> log_10(gamma)
            
            log_act_coeff=-A*sqrt(ionic_act)*(this%valence**2)/(1d0+B*this%params_act_coeff%a_TJ*sqrt(ionic_act))+&
            this%params_act_coeff%b_TJ*ionic_act
        end subroutine
    

    
end module