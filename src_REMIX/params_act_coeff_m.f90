!> This type contains the parameters that control the activity coefficient of an aqueous species
module params_act_coeff_m
    use params_aq_sol_m
    implicit none
    save
    type, public :: params_act_coeff_c !> parameters activity coefficient type
        real(kind=8) :: ion_size_param !> effective diameter hydrated ion (a_0)
        real(kind=8) :: a_TJ !> ion-specific fit parameter (Truesdell-Jones)
        real(kind=8) :: b_TJ !> ion-specific fit parameter (Truesdell-Jones)
        real(kind=8) :: b_dot !> add-on term b_dot (Helgeson)
        real(kind=8) :: alpha=0d0 !> first term
        real(kind=8) :: beta=0d0 !> second term
        real(kind=8) :: gamma=0d0 !> third term
    contains
    !> Set
        procedure, public :: set_ion_size_param
    !> Compute
        procedure, public :: compute_b_dot
        procedure, public :: compute_csts !> computes alpha and beta
    end type
    
    interface
    end interface
    
    contains
        subroutine set_ion_size_param(this,ion_size_param)
            implicit none
            class(params_act_coeff_c) :: this
            real(kind=8), intent(in) :: ion_size_param
            this%ion_size_param=ion_size_param
        end subroutine
        
        subroutine compute_b_dot(this,temp) !> Helgeson (1969)
            implicit none
            class(params_act_coeff_c) :: this
            real(kind=8), intent(in) :: temp
            this%b_dot=15698.4/temp + 41.8808*log(temp) - 0.0367626*temp - 974168/(temp**2) - 268.902
        end subroutine
        
        subroutine compute_csts(this,z,params_aq_sol,model) !> computes alpha and beta
            implicit none
            class(params_act_coeff_c) :: this
            integer(kind=4), intent(in) :: z !> valence
            class(params_aq_sol_t), intent(in) :: params_aq_sol
            integer(kind=4), intent(in) :: model
            select case (model)
            case(0) !> ideal
                this%alpha=0d0
                this%gamma=0d0
            case(1) !> Debye-H�ckel restricted
                this%alpha=1d0
                this%beta=0d0
                this%gamma=0d0
            case(2) !> Debye-H�ckel extended
                this%alpha=1d0
                this%beta=this%ion_size_param*params_aq_sol%B
                this%gamma=0d0
            case(3) !> Davies
                this%alpha=1d0
                this%beta=1d0
                this%gamma=0.3*params_aq_sol%A*(z**2)
            case(4) !> Truesdell-Jones
                this%alpha=1d0
                this%beta=params_aq_sol%B*this%a_TJ
                this%gamma=this%b_TJ
            case(5) !> b-dot
                this%alpha=1d0
                this%beta=params_aq_sol%B*this%a_TJ
                this%gamma=this%b_dot
            case default
                error stop "Model not implemented yet"
            end select
        end subroutine
end module