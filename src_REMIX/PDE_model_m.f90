module PDE_model_m
    use PDE_m
    implicit none
    save
    type, public :: PDE_model_c
        class(PDE_1D_c), pointer :: PDE
    contains
        procedure, public :: set_PDE
    end type
    
    contains
        
        subroutine set_PDE(this,PDE)
            implicit none
            class(PDE_model_c) :: this
            class(PDE_1D_c), intent(in), target :: PDE
            this%PDE=>PDE
        end subroutine
end module