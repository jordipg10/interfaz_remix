module diffusion_transient_m
    use PDE_transient_m, only: PDE_1D_transient_c, time_discr_homog_c, time_discr_heterog_c
    use diff_stab_params_m, only: stab_params_diff_c
    use diffusion_m, only: diffusion_1D_c
    use metodos_sist_lin_m!, only: tridiag_matrix_c, diag_matrix_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c, spatial_discr_c
    use conc_m, only: conc_c
    use diff_props_heterog_m, only: diff_props_heterog_c
    implicit none
    save
    type, public, extends(PDE_1D_transient_c) :: diffusion_1D_transient_c !> 1D transient diffusion subclass
        real(kind=8), allocatable :: conc_init(:)               !> initial concentration (c_0)
        type(conc_c) :: conc_obj                                !> 1D concentration class object
        type(diff_props_heterog_c) :: diff_props_heterog        !> properties
        type(stab_params_diff_c) :: stab_params_diff            !> stability parameters
    contains
        procedure, public :: allocate_conc_init
        procedure, public :: set_conc_init
        !procedure, public :: set_conc_ext
        !procedure, public :: set_conc_r_flag=>set_conc_r_flag_diff
        procedure, public :: compute_trans_mat_PDE=>compute_trans_mat_diff
        procedure, public :: compute_rech_mat_PDE=>compute_rech_mat_diff_trans
        procedure, public :: compute_F_mat_PDE=>compute_F_mat_diff
        procedure, public :: initialise_PDE=>initialise_diffusion_transient
        procedure, public :: set_stab_params_diff
        !procedure, public :: update_conc_ext
        procedure, public :: prod_total_conc
        procedure, public :: write_PDE_1D=>write_diffusion_transient_1D
        procedure, public :: set_diff_props_heterog
        procedure, public :: fund_sol_diff_eqn_1D
        procedure, public :: compute_source_term_PDE=>compute_source_term_PDE_diff_trans
        procedure, public :: solve_diff_EI_Delta_t_homog
        procedure, public :: solve_diff_RKF45
        procedure, public :: solve_diff_EE_Delta_t_heterog
        procedure, public :: solve_PDE_1D=>solve_diff_trans_1D
    end type
!****************************************************************************************************************************************************
    interface
        subroutine compute_F_mat_diff(this)
            import diffusion_1D_transient_c
            implicit none
            class(diffusion_1D_transient_c) :: this
        end subroutine
        
        subroutine compute_trans_mat_diff(this)
            import diffusion_1D_transient_c
            implicit none
            class(diffusion_1D_transient_c) :: this
        end subroutine
        
        subroutine initialise_diffusion_transient(this,root)
            import diffusion_1D_transient_c
            implicit none
            class(diffusion_1D_transient_c) :: this
            character(len=*), intent(in) :: root !> root name for output files
        end subroutine
        
        
        
        subroutine prod_total_conc(this,A_mat,time)
            import diffusion_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(diffusion_1D_transient_c) :: this
            real(kind=8), intent(in), optional :: time
            class(tridiag_matrix_c), intent(in) :: A_mat
        end subroutine
        
        subroutine write_diffusion_transient_1D(this,root,Time_out,output)
            import diffusion_1D_transient_c
            !import props_c
            implicit none
            class(diffusion_1D_transient_c), intent(in) :: this
            character(len=*), intent(in) :: root !> root name for output files
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(in) :: output(:,:)
        end subroutine
        
        subroutine solve_write_diffusion_transient(this,Time_out)
            import diffusion_1D_transient_c
            implicit none
            class(diffusion_1D_transient_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
        end subroutine
        
        subroutine solve_diff_trans_1D(this,Time_out,output)
            import diffusion_1D_transient_c
            class(diffusion_1D_transient_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        function compute_c_tilde(this,j,conc,conc_r,mixing_ratios) result(c_tilde)
            import diffusion_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(diffusion_1D_transient_c), intent(in) :: this
            integer(kind=4), intent(in) :: j
            real(kind=8), intent(in) :: conc(:,:)
            real(kind=8), intent(in) :: conc_r(:)
            class(tridiag_matrix_c), intent(in) :: mixing_ratios
            real(kind=8), allocatable :: c_tilde(:)
        end function
        
        subroutine solve_diff_EI_Delta_t_homog(this,theta,Time_out,output)
            !> Solves 1D transient PDE with homogeneous time step using Lagr explicit method 
    
            !> this: transient PDE object
            !> Time_out: output time values
            !> output: concentration vs time output
    
            !> Results at all intermediate steps are written in binary mode in file conc_binary_EE.txt
    
            !use BCs_subroutines_m, only: Dirichlet_BCs_PDE, Neumann_homog_BCs, Robin_Neumann_homog_BCs
            import diffusion_1D_transient_c
    
            !> Variables
            class(diffusion_1D_transient_c) :: this
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        subroutine solve_diff_RKF45(this,Delta_t_init,tolerance)
        import diffusion_1D_transient_c
        class(diffusion_1D_transient_c) :: this
        real(kind=8), intent(in) :: Delta_t_init
        real(kind=8), intent(in) :: tolerance
        end subroutine
        
        subroutine solve_diff_EE_Delta_t_heterog(this,Time_out,output)
        import diffusion_1D_transient_c
        class(diffusion_1D_transient_c) :: this
        real(kind=8), intent(in) :: Time_out(:)
        real(kind=8), intent(out) :: output(:,:)
        end subroutine

    end interface
!****************************************************************************************************************************************************
    contains
        subroutine set_conc_init(this,conc_init)
            implicit none
            class(diffusion_1D_transient_c) :: this
            real(kind=8), intent(in) :: conc_init(:)
            if (this%spatial_discr%Num_targets_defined.eqv..true.) then
                if (size(conc_init)/=this%spatial_discr%Num_targets) error stop "Dimension error in initial concentration"
            else
                this%spatial_discr%Num_targets=size(conc_init)
                this%spatial_discr%Num_targets_defined=.true.
            end if
            this%conc_init=conc_init
        end subroutine
    
        subroutine set_stab_params_diff(this,stab_params_diff)
            implicit none
            class(diffusion_1D_transient_c) :: this
            type(stab_params_diff_c), intent(in) :: stab_params_diff
            this%stab_params_diff=stab_params_diff
        end subroutine
        
        !subroutine set_conc_ext(this,conc_ext)
        !    class(diffusion_1D_transient_c) :: this
        !    real(kind=8), intent(in) :: conc_ext(:)
        !    if (size(conc_ext)/=this%spatial_discr%Num_targets) error stop "Dimension error in external concentration"
        !    this%conc_ext=conc_ext
        !end subroutine 
        
        !subroutine update_conc_ext(this,conc_ext_new)
        !    implicit none
        !    class(diffusion_1D_transient_c) :: this
        !    real(kind=8), intent(in) :: conc_ext_new
        !    this%conc_ext=conc_ext_new
        !end subroutine
        
        !subroutine set_conc_r_flag_diff(this)
        !    implicit none
        !    class(diffusion_1D_transient_c) :: this
        !    integer(kind=4) :: i
        !    allocate(this%conc_r_flag(this%spatial_discr%Num_targets))
        !    this%conc_r_flag=0
        !    do i=1,this%spatial_discr%Num_targets
        !        if (this%diff_props_heterog%source_term(i)>0) then
        !            this%conc_r_flag(i)=1
        !        end if
        !    end do
        !end subroutine
        !
        subroutine set_diff_props_heterog(this,diff_props_heterog)
            implicit none
            class(diffusion_1D_transient_c) :: this
            class(diff_props_heterog_c), intent(in) :: diff_props_heterog
            this%diff_props_heterog=diff_props_heterog
        end subroutine

        function fund_sol_diff_eqn_1D(this,M,Delta_x,x,mu,t) result(conc)
        !> Fundamental solution of diffusion equation in 1D
            implicit none
            class(diffusion_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: M
            real(kind=8), intent(in) :: Delta_x
            real(kind=8), intent(in) :: x
            real(kind=8), intent(in) :: mu
            real(kind=8), intent(in) :: t
            !real(kind=8), intent(in) :: phi
            !real(kind=8), intent(in) :: D
            real(kind=8) :: conc
            
            real(kind=8), parameter :: pi=4d0*atan(1d0)
            real(kind=8), parameter :: eps=1d-12
            
            !if (mod(n,2).eqv.0) then
            !    mu=(n-1)/2d0
            !else
            !    mu=floor(n/2d0)
            !end if
            
            if (abs(t)<eps) then
                if (abs(x-mu)<eps) then
                    conc=M/Delta_x
                else
                    conc=0d0
                end if
            else
                conc=(M/sqrt(4d0*pi*this%diff_props_heterog%dispersion(1)*t))*exp(-(25d-2*(x-mu)**2)/(&
                    this%diff_props_heterog%dispersion(1)*t))
            end if
        end function
        
        subroutine allocate_conc_init(this)
            implicit none
            class(diffusion_1D_transient_c) :: this
            allocate(this%conc_init(this%spatial_discr%Num_targets))
        end subroutine
        
        subroutine compute_source_term_PDE_diff_trans(this)
        implicit none
        class(diffusion_1D_transient_c) :: this
        !this%source_term_PDE=this%rech_mat%diag*this%diff%conc_ext
        this%source_term_PDE(1)=this%bd_mat(1)*this%BCs%conc_inf
        this%source_term_PDE(this%spatial_discr%Num_targets)=this%bd_mat(2)*this%BCs%conc_out
        end subroutine
        
        subroutine compute_rech_mat_diff_trans(this)
        implicit none
        class(diffusion_1D_transient_c) :: this
        !> $R$ matrix for recharge
        !> $R=diag(r)$
        !this%rech_mat%diag=this%diff%diff_props_heterog%source_term
        end subroutine
end module