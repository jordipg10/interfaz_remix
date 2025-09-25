!> Reaction module
!! This class contains the attributes of a reaction
module reaction_m
    use species_m, only: species_c
    implicit none
    save
    type, public, abstract :: reaction_c !> reaction superclass
        class(species_c), allocatable :: species(:) !> species involved
        real(kind=8), allocatable :: stoichiometry(:) !> same dimension and order as "species" attribute
        character(len=256) :: name !> name of reaction (eg. denitrification)
        integer(kind=4) :: num_species !> number of species involved
        integer(kind=4) :: react_type   !> 1: aqueous complex
                                        !! 2: mineral dissolution/precipitation
                                        !! 3: cation exchange
                                        !! 4: redox
                                        !! 5: linear
                                        !! 6: gas
                                        !! 7: nonlinear
        real(kind=8) :: eq_cst !> equilibrium constant
        real(kind=8), allocatable :: coeffs_logK_T(:) !> 6 coefficients to compute analytical expression of logK(T)
        real(kind=8) :: delta_h !> enthalpy
    contains
    !> Set
        procedure, public :: set_stoichiometry
        procedure, public :: set_all_species
        procedure, public :: set_species_names
        procedure, public :: set_single_species
        procedure, public :: set_eq_cst
        procedure, public :: set_delta_h
        procedure, public :: set_react_name
        procedure, public :: set_react_type
    !> Allocate/deallocate
        procedure, public :: allocate_reaction
        procedure, public :: deallocate_reaction
    !> Is/Are
        procedure, public :: is_species_in_react
        procedure, public :: are_species_in_react
    !> Compute
        procedure, public :: compute_logK_dep_T
    !> Write
        procedure, public :: write_reaction=>write_reaction_sup
    !> Others
        procedure, public :: copy_attributes
        procedure, public :: change_sign_stoichiometry
    end type
    
    abstract interface
       
    end interface
    
    contains
       subroutine allocate_reaction(this,num_species)
            implicit none
            class(reaction_c) :: this
            integer(kind=4), intent(in), optional :: num_species
            if (present(num_species)) then
                this%num_species=num_species
            end if
            allocate(this%stoichiometry(this%num_species),this%species(this%num_species),this%coeffs_logK_T(6))
       end subroutine
       
       subroutine deallocate_reaction(this)
            implicit none
            class(reaction_c) :: this
            deallocate(this%stoichiometry,this%species,this%coeffs_logK_T)
       end subroutine
       
       subroutine set_stoichiometry(this,stoichiometry)
            implicit none
            class(reaction_c) :: this
            real(kind=8), intent(in) :: stoichiometry(:)
            if (allocated(this%species) .and. size(stoichiometry)/=this%num_species) then
                error stop "Dimension error in 'set_stoichiometry'"
            end if
            this%stoichiometry=stoichiometry
       end subroutine
       
       subroutine set_all_species(this,species)
            implicit none
            class(reaction_c) :: this
            class(species_c), intent(in) :: species(:)
            this%species=species
       end subroutine
       
       subroutine set_single_species(this,species,index)
            implicit none
            class(reaction_c) :: this
            class(species_c), intent(in) :: species
            integer(kind=4), intent(in) :: index
            if (index<1 .or. index>this%num_species) error stop
            call this%species(index)%assign_species(species)
       end subroutine

       subroutine set_species_names(this,species_names)
            implicit none
            class(reaction_c) :: this
            character(len=*), intent(in) :: species_names(:)
            integer(kind=4) :: i
            if (size(species_names).ne.this%num_species) error stop "Dimension error in 'set_species_names'"
            do i=1,this%num_species
                call this%species(i)%set_name(species_names(i))
            end do
       end subroutine
              
       subroutine set_eq_cst(this,eq_cst)
            implicit none
            class(reaction_c) :: this
            real(kind=8), intent(in) :: eq_cst
            this%eq_cst=eq_cst
       end subroutine
       
       subroutine set_delta_h(this,delta_h)
            implicit none
            class(reaction_c) :: this
            real(kind=8), intent(in) :: delta_h
            this%delta_h=delta_h
       end subroutine
       
       subroutine set_react_name(this,name)
            implicit none
            class(reaction_c) :: this
            character(len=*), intent(in) :: name
            this%name=name
       end subroutine
       
       subroutine set_react_type(this,react_type)
            implicit none
            class(reaction_c) :: this
            integer(kind=4), intent(in) :: react_type
            this%react_type=react_type
       end subroutine
       
       subroutine is_species_in_react(this,species,flag,species_ind) !> checks if a species is involved in reaction
            implicit none
            class(reaction_c), intent(in) :: this
            class(species_c), intent(in) :: species
            logical, intent(out) :: flag
            integer(kind=4), intent(out), optional :: species_ind !> species index in reaction
            
            integer(kind=4) :: i
            
            flag=.false.
            if (present(species_ind)) then
                species_ind=0
            end if
            do i=1,this%num_species
                if (species%name==this%species(i)%name) then
                    flag=.true.
                    if (present(species_ind)) then
                        species_ind=i
                    end if
                    exit
                end if
            end do
       end subroutine
       
       subroutine are_species_in_react(this,species,flag,species_ind) !> checks if a set of species participates in a reaction
            implicit none
            class(reaction_c), intent(in) :: this
            class(species_c), intent(in) :: species(:)
            logical, intent(out) :: flag
            integer(kind=4), intent(out) :: species_ind(:) !> species indices in reaction (already allocated)
            
            integer(kind=4) :: i,j
            
            if (size(species)/=size(species_ind)) then
                error stop "Dimension error in subroutine 'are_species_in_react'"
            end if
            flag=.true.
            i=1
            j=1
            do
                if (species(j)%name==this%species(i)%name) then
                    species_ind(j)=i
                    if (j<size(species)) then
                        j=j+1
                        i=1
                    else
                        exit
                    end if
                else if (i<this%num_species) then
                    i=i+1
                else
                    flag=.false.
                    exit
                end if
            end do
       end subroutine
       
        subroutine copy_attributes(this,reaction)
            implicit none
            class(reaction_c) :: this
            class(reaction_c), intent(in) :: reaction
            this%species=reaction%species
            this%stoichiometry=reaction%stoichiometry
            this%name=reaction%name
            this%num_species=reaction%num_species
            this%react_type=reaction%react_type
            this%delta_h=reaction%delta_h
            this%eq_cst=reaction%eq_cst
            this%coeffs_logK_T=reaction%coeffs_logK_T
        end subroutine
        
        subroutine compute_logK_dep_T(this,temp) !> based on PHREEQC
            implicit none
            class(reaction_c) :: this
            real(kind=8), intent(in) :: temp !> temperature of solution
            
            real(kind=8) :: logK
            
            if (size(this%coeffs_logK_T)/=6) error stop "There must be 6 coefficients to compute logK(T)"
            logK=this%coeffs_logK_T(1)+this%coeffs_logK_T(2)*temp+this%coeffs_logK_T(3)/temp+this%coeffs_logK_T(4)*log10(temp)+&
            this%coeffs_logK_T(5)/(temp**2)+this%coeffs_logK_T(6)*temp**2
            this%eq_cst=10**logK
        end subroutine
        
        subroutine allocate_coeffs_logK_T(this) !> based on PHREEQC
            implicit none
            class(reaction_c) :: this 
            allocate(this%coeffs_logK_T(6))
        end subroutine
        
        subroutine change_sign_stoichiometry(this)
            implicit none
            class(reaction_c) :: this  
            integer(kind=4) :: i
            do i=1,this%num_species
                this%stoichiometry(i)=-this%stoichiometry(i)
            end do
        end subroutine
        
        subroutine write_reaction_sup(this,unit)
            implicit none
            class(reaction_c) :: this  
            integer(kind=4), intent(in) :: unit
            write(unit,"(10x,A30/)") this%name
            write(unit,"(10x,I5/)") this%react_type
            write(unit,"(10x,I5/)") this%num_species
            write(unit,"(10x,ES15.5/)") this%eq_cst
       end subroutine
end module