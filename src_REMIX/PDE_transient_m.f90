!> Transient PDE module
!> $F dc/dt=Tc+Ec_ext+Bc_bd$, where $c$ is the state variable
!> $\theta$-method: A * c^{k+1}= X*c^k + Y*c_rech + Z*c_bd
module PDE_transient_m
    use PDE_m
    use time_discr_m, only: time_discr_c, time_discr_homog_c, time_discr_heterog_c
    use char_params_m, only: char_params_c
    implicit none
    save
    type, public, abstract, extends(PDE_1D_c) :: PDE_1D_transient_c !> 1D transient PDE subclass
        class(time_discr_c), pointer :: time_discr !> time discretisation (polymorphic variable)
        !class(char_params_c), pointer :: char_params !> characteristic parameters (polymorphic variable)
        type(diag_matrix_c) :: F_mat !> storage matrix
        type(tridiag_matrix_c) :: X_mat !> mixing ratios domain previous time step
        type(diag_matrix_c) :: Y_mat !> mixing ratios external
        real(kind=8), allocatable :: Z_mat(:) !> mixing ratios boundary (dim=2)
        type(tridiag_matrix_c) :: A_mat !>  mixing ratios domain next time step
        real(kind=8), allocatable :: A_mat_inv(:,:) !> inverse of A_mat
        real(kind=8), allocatable :: f_vec(:) !> independent term in linear system (f_vec=Y*c_rech + Z*c_bd)
        type(real_array_c) :: mixing_ratios_conc !> mixing ratios concentrations WMA (vale para cualquier dimension)
        type(real_array_c) :: mixing_ratios_Rk !> mixing ratios reaction amounts WMA (vale para cualquier dimension)
        type(real_array_c) :: mixing_ratios_Rk_init !> initial mixing ratios reaction amounts WMA (vale para cualquier dimension)
        type(int_array_c) :: mixing_waters_indices !> matrix that contains indices of target waters that mix with each target water
        type(int_array_c) :: mixing_waters_indices_dom !> matrix that contains indices of domain target waters that mix with each domain target water
        real(kind=8), allocatable :: mixing_ratios_mat_conc_dom(:,:) !> mixing ratios domain matrix concentrations WMA (Euler implicito)
        !real(kind=8), allocatable :: mixing_ratios_mat_Rk(:,:) !> mixing ratios matrix reaction amounts WMA (Euler implicito)
        !real(kind=8), allocatable :: mixing_ratios_mat_Rk_init(:,:) !> mixing ratios matrix reaction amounts WMA (Euler implicito)
    contains
    !> set
        procedure, public :: set_time_discr
        !procedure, public :: set_char_params
    !> Allocate
        procedure, public :: allocate_arrays_PDE_1D=>allocate_arrays_PDE_1D_trans
        procedure, public :: allocate_F_mat
        procedure, public :: allocate_X_mat
        procedure, public :: allocate_Y_mat
        procedure, public :: allocate_Z_mat
        procedure, public :: allocate_A_mat
        procedure, public :: allocate_A_mat_inv
        procedure, public :: allocate_f_vec
        procedure, public :: allocate_mixing_ratios
        procedure, public :: allocate_mixing_ratios_mat_conc_dom
        procedure, public :: allocate_mixing_waters_indices
    !> Computations
        procedure(compute_F_mat_PDE), public, deferred :: compute_F_mat_PDE
        procedure, public :: compute_E_mat
        procedure, public :: compute_X_mat
        procedure, public :: compute_Y_mat
        procedure, public :: compute_Z_mat
        procedure, public :: compute_A_mat
        procedure, public :: compute_lumped_A_mat
        !procedure, public :: compute_mixing_ratios_Delta_t_homog
        procedure, public :: compute_f_vec
        procedure, public :: compute_b_vec_lin_syst
        procedure, public :: compute_A_mat_ODE
        procedure, public :: compute_b_ODE     
        !procedure, public :: solve_PDE_EE_Delta_t_homog
        !procedure, public :: solve_PDE_EE_Delta_t_heterog
        !procedure, public :: solve_PDE_EI_Delta_t_homog
        !procedure, public :: solve_PDE_RKF45
        procedure, public :: compute_k_RKF45
    end type
!*****************************************************************************************************************************
    abstract interface
        
        subroutine compute_F_mat_PDE(this)
            import PDE_1D_transient_c
            implicit none
            class(PDE_1D_transient_c) :: this
        end subroutine
        
    end interface
!*****************************************************************************************************************************
    interface
        
        subroutine compute_E_mat(this,E_mat,k)
            import PDE_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c) :: this
            type(tridiag_matrix_c), intent(out) :: E_mat
            integer(kind=4), intent(in), optional :: k
        end subroutine
        
        
        
        subroutine compute_A_mat(this,theta,E_mat)
            import PDE_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: theta
            class(tridiag_matrix_c), intent(in) :: E_mat
        end subroutine
        
        subroutine compute_lumped_A_mat(this,A_mat_lumped)
            import PDE_1D_transient_c
            import diag_matrix_c
            implicit none
            class(PDE_1D_transient_c), intent(in) :: this
            type(diag_matrix_c), intent(out) :: A_mat_lumped !> must be allocated
        end subroutine
        
        !subroutine compute_mixing_ratios_Delta_t_homog(this,A_mat_lumped)
        !    import PDE_1D_transient_c
        !    import diag_matrix_c
        !    implicit none
        !    class(PDE_1D_transient_c) :: this
        !    !real(kind=8), intent(in) :: theta
        !    type(diag_matrix_c), intent(out), optional :: A_mat_lumped
        !end subroutine
        
        subroutine compute_f_vec(this,k)
            import PDE_1D_transient_c
            implicit none
            class(PDE_1D_transient_c) :: this
            integer(kind=4), intent(in), optional :: k
        end subroutine
        
        subroutine compute_b_vec_lin_syst(this,theta,conc_old,b,k)
            import PDE_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(inout) :: b(:)
            integer(kind=4), intent(in), optional :: k
        end subroutine
        
        subroutine compute_A_mat_ODE(this,A_mat)
            import PDE_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c), intent(in) :: this
            type(tridiag_matrix_c), intent(out) :: A_mat
        end subroutine
        
        
        function compute_b_ODE(this) result(b)
            import PDE_1D_transient_c
            implicit none
            class(PDE_1D_transient_c), intent(in) :: this
            real(kind=8), allocatable :: b(:)
        end function
        
        subroutine solve_PDE_EE_Delta_t_homog(this,Time_out,output)
            import PDE_1D_transient_c
            class(PDE_1D_transient_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        subroutine solve_PDE_EE_Delta_t_heterog(this,Time_out,output)
            import PDE_1D_transient_c
            class(PDE_1D_transient_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        subroutine solve_PDE_EI_Delta_t_homog(this,theta,Time_out,output)
            import PDE_1D_transient_c
            class(PDE_1D_transient_c) :: this
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        subroutine solve_PDE_RKF45(this,Delta_t_init,tolerance)
            import PDE_1D_transient_c
            class(PDE_1D_transient_c) :: this
            real(kind=8), intent(in) :: Delta_t_init
            real(kind=8), intent(in) :: tolerance
        end subroutine
        
        subroutine update_time_step_RKF45(Delta_t_old,tolerance,conc_RK4,conc_RK5,Delta_t_new)
            implicit none
            real(kind=8), intent(in) :: Delta_t_old
            real(kind=8), intent(in) :: tolerance
            real(kind=8), intent(in) :: conc_RK4(:)
            real(kind=8), intent(in) :: conc_RK5(:)
            real(kind=8), intent(out) :: Delta_t_new
        end subroutine
        
        function compute_k_RKF45(this,Delta_t,conc_RK4) result(k)
            import PDE_1D_transient_c
            implicit none
            class(PDE_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: conc_RK4(:)
            real(kind=8), allocatable :: k(:,:)
        end function
    end interface
!*****************************************************************************************************************************
    contains
        subroutine set_time_discr(this,time_discr_obj)
            implicit none
            class(PDE_1D_transient_c) :: this
            class(time_discr_c), intent(in), target :: time_discr_obj
            this%time_discr=>time_discr_obj
        end subroutine
        
        !subroutine set_char_params(this,char_params_obj)
        !    implicit none
        !    class(PDE_1D_transient_c) :: this
        !    class(char_params_c), intent(in), target :: char_params_obj
        !    this%char_params=>char_params_obj
        !end subroutine
        
        subroutine allocate_F_mat(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            !print *, "Allocating F_mat of size ", this%spatial_discr%Num_targets
            call this%F_mat%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
        subroutine allocate_X_mat(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            call this%X_mat%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
        subroutine allocate_Y_mat(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            call this%Y_mat%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
        subroutine allocate_Z_mat(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            allocate(this%Z_mat(2))
        end subroutine
        
        subroutine allocate_A_mat(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            call this%A_mat%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
        subroutine allocate_f_vec(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            allocate(this%f_vec(this%spatial_discr%Num_targets))
        end subroutine
        
        subroutine allocate_mixing_ratios(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            call this%mixing_ratios_conc%allocate_array(this%spatial_discr%Num_targets)
            call this%mixing_ratios_Rk%allocate_array(this%spatial_discr%Num_targets)
            call this%mixing_ratios_Rk_init%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
        subroutine allocate_mixing_ratios_mat_conc_dom(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            allocate(this%mixing_ratios_mat_conc_dom(this%spatial_discr%Num_targets,this%spatial_discr%Num_targets))
            !call this%mixing_ratios_Rk%allocate_array(this%spatial_discr%Num_targets)
            !call this%mixing_ratios_Rk_init%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
         subroutine allocate_A_mat_inv(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            allocate(this%A_mat_inv(this%spatial_discr%Num_targets,this%spatial_discr%Num_targets))
            !call this%mixing_ratios_conc%allocate_array(this%spatial_discr%Num_targets)
            !call this%mixing_ratios_Rk%allocate_array(this%spatial_discr%Num_targets)
            !call this%mixing_ratios_Rk_init%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
        subroutine allocate_arrays_PDE_1D_trans(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            call allocate_arrays_PDE_1D_stat(this)
            call this%allocate_F_mat()
            call this%allocate_X_mat()
            call this%allocate_Y_mat()
            call this%allocate_Z_mat()
            call this%allocate_A_mat()
            call this%allocate_f_vec()
        end subroutine

        
        subroutine allocate_mixing_waters_indices(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            integer(kind=4) :: i
            call this%mixing_waters_indices%allocate_array(this%mixing_ratios_conc%num_cols)
            call this%mixing_waters_indices_dom%allocate_array(this%mixing_ratios_conc%num_cols)
            do i=1,this%mixing_waters_indices%num_cols
                call this%mixing_waters_indices%cols(i)%allocate_vector(this%mixing_ratios_conc%cols(i)%dim+2)
                call this%mixing_waters_indices_dom%cols(i)%allocate_vector(this%spatial_discr%Num_targets+2)
            end do
        end subroutine

        subroutine compute_X_mat(this,theta,E_mat)
            !import PDE_1D_transient_c
            !import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c) :: this
            real(kind=8), intent(in) :: theta
            class(tridiag_matrix_c), intent(in) :: E_mat
            
            integer(kind=4) :: n
            real(kind=8) :: B_norm_inf,B_norm_1
        
            !call this%X_mat%allocate_array(n)
            this%X_mat%sub=(1d0-theta)*E_mat%sub
            this%X_mat%diag=1d0+(1d0-theta)*E_mat%diag
            this%X_mat%super=(1d0-theta)*E_mat%super
        end subroutine
        
        subroutine compute_Y_mat(this,k)
            !import PDE_1D_transient_c
            !import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c) :: this
            !real(kind=8), intent(in) :: theta
            integer(kind=4), intent(in), optional :: k
            !class(tridiag_matrix_c), intent(in) :: E_mat
            
            integer(kind=4) :: n
            real(kind=8) :: B_norm_inf,B_norm_1
        
            select type (time_discr=>this%time_discr)
            type is (time_discr_homog_c)
                if (this%dimless) then
                    this%Y_mat%diag=this%rech_mat%diag*time_discr%Delta_t_D
                else
                    this%Y_mat%diag=this%rech_mat%diag*time_discr%Delta_t
                end if
                !E_mat%super=E_mat%super*time_discr%Delta_t
            type is (time_discr_heterog_c)
                if (this%dimless) then
                    this%Y_mat%diag=this%rech_mat%diag*time_discr%Delta_t_D(k)
                else
                    this%Y_mat%diag=this%rech_mat%diag*time_discr%Delta_t(k)
                end if
            end select
            this%Y_mat%diag=this%Y_mat%diag/this%F_mat%diag
            !call this%X_mat%allocate_array(n)
            !this%Y_mat%sub=(1d0-theta)*this%rech_mat%sub
            !this%Y_mat%diag=1d0+(1d0-theta)*this%rech_mat%diag
            !this%Y_mat%super=(1d0-theta)*E_mat%super
        end subroutine
       
        subroutine compute_Z_mat(this,k)
            !import PDE_1D_transient_c
            !import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c) :: this
            !real(kind=8), intent(in) :: theta
            integer(kind=4), intent(in), optional :: k
            !class(tridiag_matrix_c), intent(in) :: E_mat
            
            integer(kind=4) :: n
            real(kind=8) :: B_norm_inf,B_norm_1
        
            select type (time_discr=>this%time_discr)
            type is (time_discr_homog_c)
                if (this%dimless) then
                    this%Z_mat=this%bd_mat*time_discr%Delta_t_D
                else
                    this%Z_mat=this%bd_mat*time_discr%Delta_t
                end if
                !this%Z_mat=this%bd_mat*time_discr%Delta_t
            type is (time_discr_heterog_c)
                if (this%dimless) then
                    this%Z_mat=this%bd_mat*time_discr%Delta_t_D(k)
                else
                    this%Z_mat=this%bd_mat*time_discr%Delta_t(k)
                end if
            end select
            this%Z_mat(1)=this%Z_mat(1)/this%F_mat%diag(1)
            this%Z_mat(2)=this%Z_mat(2)/this%F_mat%diag(this%spatial_discr%num_targets)
            !call this%X_mat%allocate_array(n)
            !this%Y_mat%sub=(1d0-theta)*this%rech_mat%sub
            !this%Y_mat%diag=1d0+(1d0-theta)*this%rech_mat%diag
            !this%Y_mat%super=(1d0-theta)*E_mat%super
        end subroutine

end module