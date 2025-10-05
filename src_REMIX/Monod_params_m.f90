!> This subclass contains the parameters of a Monod reaction rate
module Monod_params_m
    use kin_params_m, only: kin_params_c
    implicit none
    save
    type, public, extends(kin_params_c) :: Monod_params_c
        integer(kind=4) :: num_terms !> number of terms
        integer(kind=4) :: n_inh !> number of inhibitors
        real(kind=8), allocatable :: conc_thr(:) !> threshold concentration electron donors & acceptors
        real(kind=8), allocatable :: k_M(:) !> half-saturation constants
        real(kind=8), allocatable :: k_inh(:) !> inhibition constants
        real(kind=8), allocatable :: conc_thr_inh(:) !> threshold concentration inhibitors
        
    contains
        procedure, public :: allocate_k_inh
        procedure, public :: allocate_k_M
        procedure, public :: compute_num_terms
    end type
        
    interface
     
    end interface
    
    contains
      
        subroutine allocate_k_M(this)
            implicit none
            class(Monod_params_c) :: this
            allocate(this%k_M(2))
        end subroutine
        
        subroutine allocate_k_inh(this,n_inh)
            implicit none
            class(Monod_params_c) :: this
            integer(kind=4), intent(in), optional :: n_inh
            if (present(n_inh)) then
                if (n_inh<0) then
                    error stop "Number of inhibitors must be positive"
                else
                    this%n_inh=n_inh
                end if
            end if
            allocate(this%k_inh(this%n_inh))
        end subroutine
        
        subroutine compute_num_terms(this)
            implicit none
            class(Monod_params_c) :: this
            this%num_terms=this%n_inh+2
        end subroutine
        
end module