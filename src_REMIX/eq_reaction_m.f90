!> Equilibrium reaction subclass
!! This class contains attributes of equilibrium reactions
module eq_reaction_m
    use reaction_m
    use aq_species_m
    use phase_m
    implicit none
    save
    type, public, extends(reaction_c) :: eq_reaction_c !> equilibrium reaction subclass
    contains
        procedure, public :: read_eq_reaction
        procedure, public :: read_dissolution_react_PHREEQC
        procedure, public :: read_association_react_PHREEQC
        procedure, public :: read_exchange_react_PHREEQC
    end type
    
     interface
       
        subroutine read_eq_reaction(this,species,filename)
            import eq_reaction_c
            import species_c
            implicit none
            class(eq_reaction_c) :: this
            class(species_c), intent(in) :: species
            character(len=*), intent(in) :: filename
        end subroutine
        
        subroutine read_dissolution_react_PHREEQC(this,string,phase)
            import eq_reaction_c
            import phase_c
            implicit none
            class(eq_reaction_c) :: this
            character(len=*), intent(in) :: string
            class(phase_c), intent(inout), optional :: phase !> defined phase
        end subroutine
        
        subroutine read_association_react_PHREEQC(this,string,prim_flag,defined_species)
            import eq_reaction_c
            import aq_species_c
            implicit none
            class(eq_reaction_c) :: this
            character(len=*), intent(in) :: string !> association reaction
            logical, intent(out) :: prim_flag !> TRUE if species is primary
            type(aq_species_c), intent(out), optional :: defined_species
        end subroutine
        
        subroutine read_exchange_react_PHREEQC(this,string,prim_flag,defined_species)
            import eq_reaction_c
            import species_c
            implicit none
            class(eq_reaction_c) :: this
            character(len=*), intent(in) :: string !> half reaction
            logical, intent(out) :: prim_flag !> TRUE if species is primary
            type(species_c), intent(out), optional :: defined_species
        end subroutine
        
        !  subroutine react_rate_bin_syst_eq_1D(this,u,du_dx,D,phi,r_eq)
        !     import eq_reaction_c
        !     implicit none
        !     class(eq_reaction_c) :: this
        !     real(kind=8), intent(in) :: u
        !     real(kind=8), intent(in) :: du_dx
        !     real(kind=8), intent(in) :: D
        !     real(kind=8), intent(in) :: phi
        !     real(kind=8), intent(out) :: r_eq
        !  end subroutine
    end interface
    
    contains
       
    
    
    
    
    subroutine append_eq_reaction(this,eq_reactions) !> appends eq
        implicit none
        class(eq_reaction_c), intent(in) :: this
        type(eq_reaction_c), intent(inout), allocatable :: eq_reactions(:) !> must be already allocated
        
        type(eq_reaction_c), allocatable :: aux(:)
        
        aux=eq_reactions
        deallocate(eq_reactions)
        if (size(aux)==0) then
            allocate(eq_reactions(1))
            eq_reactions(1)=this
        else
            allocate(eq_reactions(size(aux)+1))
            eq_reactions(1:size(aux))=aux
            eq_reactions(size(eq_reactions))=this
        end if
    end subroutine

    !> This subroutine computes the reaction rate of a binary system in equilibrium in 1D (De Simoni et al, 2005)
subroutine react_rate_bin_syst_eq_1D(this,u,du_dx,D,phi,r_eq)
    implicit none
    class(eq_reaction_c) :: this
    real(kind=8), intent(in) :: u !> concentration component at this target
    real(kind=8), intent(in) :: du_dx !> gradient of u at this target
    real(kind=8), intent(in) :: D !> dispersion coefficient
    real(kind=8), intent(in) :: phi !> porosity
    real(kind=8), intent(out) :: r_eq !> reaction rate
    
    real(kind=8) :: d2c2nc_du2
    
    d2c2nc_du2=2*this%eq_cst/((u**2+4*this%eq_cst)**1.5) !> second derivative of secondary concentrations with respect to component concentrations
    
    r_eq=phi*d2c2nc_du2*D*du_dx**2
   
end subroutine
end module