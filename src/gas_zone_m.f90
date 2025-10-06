module gas_zone_m !> creo que no es necesario
    use gas_m
    implicit none
    save
    type, public :: gas_zone_c !> gas zone class
        integer(kind=4) :: num_species
        integer(kind=4) :: num_species_eq
        type(gas_c), allocatable :: gases(:)
        type(gas_c), allocatable :: gases_eq(:)
    contains
        procedure, public :: set_num_species
        procedure, public :: set_num_species_eq
        procedure, public :: set_gases
        procedure, public :: set_gases_eq
        procedure, public :: allocate_gases
        procedure, public :: allocate_gases_eq
        
    end type
    
    interface
        
    end interface
    
    contains
    
        subroutine set_num_species(this,num_species)
            implicit none
            class(gas_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_species
            if (present(num_species)) then
                this%num_species=num_species
            else
                this%num_species=size(this%gases)
            end if
        end subroutine
        
        subroutine set_num_species_eq(this,num_species_eq)
            implicit none
            class(gas_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_species_eq
            if (present(num_species_eq)) then
                this%num_species_eq=num_species_eq
            else
                this%num_species_eq=size(this%gases_eq)
            end if
        end subroutine
        
        subroutine allocate_gases(this,num_species)
            implicit none
            class(gas_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_species
            if (present(num_species)) then
                this%num_species=num_species
            end if
            allocate(this%gases(this%num_species))
        end subroutine
        
        subroutine allocate_gases_eq(this,num_species_eq)
            implicit none
            class(gas_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_species_eq
            if (present(num_species_eq)) then
                this%num_species_eq=num_species_eq
            end if
            allocate(this%gases_eq(this%num_species_eq))
        end subroutine
        
        subroutine set_gases(this,gases)
            implicit none
            class(gas_zone_c) :: this
            class(gas_c), intent(in) :: gases(:)
            this%gases=gases
        end subroutine
        
        subroutine set_gases_eq(this,gases_eq)
            implicit none
            class(gas_zone_c) :: this
            class(gas_c), intent(in) :: gases_eq(:)
            !print *, allocated(this%species)
            if (size(gases_eq)>this%num_species) then
                error stop "Number of gases in equilibrium cannot be greater than number of gases"
            else
                this%gases_eq=gases_eq
            end if
        end subroutine
end module