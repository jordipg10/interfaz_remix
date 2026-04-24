!> \file diffusion_transient_m.f90
!> \brief Transient diffusion module for 1D and 2D problems.
!> \details
!> Defines the `diffusion_1D_transient_c` and `diffusion_2D_transient_c`
!> classes extending the PDE transient base classes for solving the
!> transient diffusion equation:
!> \f[
!>   \phi \frac{\partial c}{\partial t} = \nabla \cdot (D \nabla c) + r
!> \f]
!> Supports Euler implicit, explicit, and Runge-Kutta-Fehlberg (RKF45)
!> time integration schemes. Includes computation of transfer matrices,
!> stability parameters, and mixing ratios.
!>
!> \see diffusion_m, PDE_transient_m, diff_props_heterog_m
!> \author Jordi
!> \date Unknown
!> \ingroup transport
module diffusion_transient_m
    use PDE_transient_m, only: PDE_1D_transient_c, PDE_2D_transient_c
    use time_discr_m, only: time_discr_homog_c, time_discr_heterog_c
    use diff_stab_params_m, only: stab_params_diff_c
    use diffusion_m, only: diffusion_1D_c
    use arrays_m, only: tridiag_matrix_c, diag_matrix_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    use conc_m, only: conc_c
    use diff_props_heterog_m, only: diff_props_heterog_1D_c, diff_props_heterog_2D_c
    implicit none
    save
    private
    public :: fund_sol_diff_eqn_1D
    !> \brief 1D transient diffusion class.
    !> \details Extends PDE_1D_transient_c with diffusion-specific properties,
    !> stability parameters, and solver methods (EI, EE, RKF45).
    type, public, extends(PDE_1D_transient_c) :: diffusion_1D_transient_c
        real(kind=8), allocatable :: conc_init(:)               !< [M/L^3] Initial concentration
        type(conc_c) :: conc_obj                                !< Concentration storage object
        type(diff_props_heterog_1D_c) :: diff_props_heterog     !< 1D diffusion properties (porosity, D)
        type(stab_params_diff_c) :: stab_params_diff            !< Stability parameters for time-stepping
        !> \var Lagr_flag
        !> Method selector: .true. = Lagrangian, .false. = Eulerian
        
    contains
        procedure :: allocate_conc_init
        procedure :: set_conc_init
        !> Set Lagrangian/Eulerian method flag
        
        !procedure :: set_conc_ext
        !procedure :: set_conc_r_flag=>set_conc_r_flag_diff
        procedure :: compute_trans_mat_PDE=>compute_trans_mat_diff_trans_1D
        procedure :: compute_rech_mat_PDE=>compute_rech_mat_diff_trans_1D
        procedure :: compute_F_mat_PDE=>compute_F_mat_diff_trans_1D
        procedure :: initialise_PDE=>initialise_diffusion_transient
        procedure :: set_stab_params_diff
        procedure :: compute_mixing_ratios_Delta_t_homog=>compute_mixing_ratios_Delta_t_homog_1D_diff
        procedure :: prod_total_conc
        procedure :: write_PDE=>write_diffusion_transient_1D
        procedure :: set_diff_props_heterog
        procedure :: fund_sol_diff_eqn_1D
        procedure :: compute_source_term_PDE=>compute_source_term_PDE_diff_trans
        procedure :: solve_diff_EI_Delta_t_homog
        procedure :: solve_diff_RKF45
        procedure :: solve_diff_EE_Delta_t_heterog
        procedure :: solve_PDE=>solve_diff_trans_1D
    end type

    type, public, extends(PDE_2D_transient_c) :: diffusion_2D_transient_c
        !> Array of 1D transient diffusion objects for multiple species
        type(diff_props_heterog_2D_c) :: diff_props_heterog        !> properties
        real(kind=8), allocatable :: conc_init(:,:)               !> initial concentration (c_0)
    contains
        !> Set initial concentration
        !procedure :: set_conc_init_2D
        !> Initialise 2D transient diffusion problem
        !procedure :: initialise_diffusion_transient_2D
        !> Solve 2D transient diffusion problem
        procedure :: compute_F_mat_PDE=>compute_F_mat_diff_trans_2D
        procedure :: compute_source_term_pde=>compute_source_term_diff_trans_2D
        procedure :: compute_trans_mat_PDE=>compute_trans_mat_diff_trans_2D
        procedure :: initialise_PDE=>initialise_diffusion_transient_2D
        procedure :: solve_PDE=>solve_diffusion_transient_2D
        procedure :: write_PDE=>write_diffusion_transient_2D
        procedure :: compute_rech_mat_PDE=>compute_rech_mat_diff_trans_2D
        procedure :: compute_mixing_ratios_Delta_t_homog=>compute_mixing_ratios_Delta_t_homog_2D_diff
    end type
!****************************************************************************************************************************************************
    interface
        ! subroutine compute_F_mat_diff(this)
        !     import diffusion_1D_transient_c
        !     implicit none
        !     class(diffusion_1D_transient_c) :: this
        ! end subroutine
        
        subroutine compute_trans_mat_diff_trans_1D(this)
            import diffusion_1D_transient_c
            implicit none
            class(diffusion_1D_transient_c) :: this
        end subroutine
        
        subroutine initialise_diffusion_transient(this,path,root,mesh_type)
            import diffusion_1D_transient_c
            implicit none
            class(diffusion_1D_transient_c) :: this
            character(len=*), intent(in) :: path !> Input file path
            character(len=*), intent(in) :: root !> Root name for input files
            integer(kind=4), intent(in) :: mesh_type
        end subroutine
        
        
        
        subroutine prod_total_conc(this,A_mat,time)
            import diffusion_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(diffusion_1D_transient_c) :: this
            real(kind=8), intent(in), optional :: time
            class(tridiag_matrix_c), intent(in) :: A_mat
        end subroutine
        
        subroutine write_diffusion_transient_1D(this)
            import diffusion_1D_transient_c
            !import props_c
            implicit none
            class(diffusion_1D_transient_c), intent(in) :: this
            ! character(len=*), intent(in) :: root !> root name for output files
            ! real(kind=8), intent(in) :: Time_out(:)
            ! real(kind=8), intent(in) :: output(:,:)
        end subroutine
        
        subroutine solve_write_diffusion_transient(this,Time_out)
            import diffusion_1D_transient_c
            implicit none
            class(diffusion_1D_transient_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
        end subroutine
        
        subroutine solve_diff_trans_1D(this)
            import diffusion_1D_transient_c
            class(diffusion_1D_transient_c) :: this
            !real(kind=8), intent(in) :: Time_out(:)
            !real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        function compute_c_mix(this,j,conc,conc_r,mixing_ratios) result(c_mix)
            import diffusion_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(diffusion_1D_transient_c), intent(in) :: this
            integer(kind=4), intent(in) :: j
            real(kind=8), intent(in) :: conc(:,:)
            real(kind=8), intent(in) :: conc_r(:)
            class(tridiag_matrix_c), intent(in) :: mixing_ratios
            real(kind=8), allocatable :: c_mix(:)
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
            class(diff_props_heterog_1D_c), intent(in) :: diff_props_heterog
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
                conc=(M/sqrt(4d0*pi*this%diff_props_heterog%diff_cent(1)*t))*exp(-(25d-2*(x-mu)**2)/(&
                    this%diff_props_heterog%diff_cent(1)*t))
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
        
        subroutine compute_rech_mat_diff_trans_1D(this)
        implicit none
        class(diffusion_1D_transient_c) :: this
        !> $R$ matrix for recharge
        !> $R=diag(r)$
        !this%rech_mat%diag=this%diff%diff_props_heterog%source_term
        end subroutine

        

        subroutine compute_F_mat_diff_trans_1D(this) !> diagonal matrix
        implicit none
        class(diffusion_1D_transient_c) :: this
        
        integer(kind=4) :: i,n
        real(kind=8) :: r_i
        
        n=this%spatial_discr%Num_targets
        if (this%spatial_discr%adapt_ref==1) then
            deallocate(this%F_mat%diag)
            call this%F_mat%allocate_array(n)
        end if
        this%F_mat%diag=this%diff_props_heterog%porosity
        select type (mesh=>this%spatial_discr)
        type is (spatial_discr_rad_c)
            if (mesh%dim>1) then
                !this%F_mat%diag=this%diff_props_heterog%porosity
                forall (i=1:n)
                    this%F_mat%diag(i)=this%F_mat%diag(i)*(i-5d-1)**(mesh%dim-1)
                end forall
            end if
        end select
        this%F_mat_prev%diag=this%F_mat%diag
    end subroutine

    subroutine compute_F_mat_diff_trans_2D(this) !> diagonal matrix
    implicit none   
    class(diffusion_2D_transient_c) :: this
    
    integer(kind=4) :: i,n
    n=this%spatial_discr%Num_targets
    if (this%spatial_discr%adapt_ref==1) then
        deallocate(this%F_mat%diag)
        call this%F_mat%allocate_array(n)
    end if
    this%F_mat%diag=this%diff_props_heterog%porosity
    this%F_mat_prev%diag=this%F_mat%diag
    end subroutine
    
    subroutine compute_source_term_diff_trans_2D(this)
    implicit none
    class(diffusion_2D_transient_c) :: this
    ! Implementation needed
    end subroutine
    
    subroutine compute_trans_mat_diff_trans_2D(this)
    implicit none
    class(diffusion_2D_transient_c) :: this
    ! Implementation needed
    end subroutine
    
    subroutine initialise_diffusion_transient_2D(this,path,root,mesh_type)
    implicit none
    class(diffusion_2D_transient_c) :: this
    character(len=*), intent(in) :: path !> Input file path
    character(len=*), intent(in) :: root !> Root name for input files
    integer(kind=4), intent(in) :: mesh_type
    ! Implementation needed
    end subroutine
    
    subroutine write_diffusion_transient_2D(this)
    implicit none
    class(diffusion_2D_transient_c), intent(in)  :: this
    ! Implementation needed
    end subroutine
    
    subroutine solve_diffusion_transient_2D(this)
    implicit none
    class(diffusion_2D_transient_c) :: this
    ! Implementation needed
    end subroutine
    
    subroutine compute_rech_mat_diff_trans_2D(this)
    implicit none
    class(diffusion_2D_transient_c) :: this
    ! Implementation needed
    end subroutine

    subroutine compute_mixing_ratios_Delta_t_homog_1D_diff(this)
    implicit none
    class(diffusion_1D_transient_c) :: this
    !real(kind=8), intent(in) :: theta
    !type(diag_matrix_c), intent(out), optional :: A_mat_lumped
    ! For diffusion, mixing ratios are not needed, so this can be left empty or used to compute a lumped matrix if needed for stability analysis
    end subroutine

    subroutine compute_mixing_ratios_Delta_t_homog_2D_diff(this)
    implicit none
    class(diffusion_2D_transient_c) :: this
    !real(kind=8), intent(in) :: theta
    !type(diag_matrix_c), intent(out), optional :: A_mat_lumped
    ! For diffusion, mixing ratios are not needed, so this can be left empty or used to compute a lumped matrix if needed for stability analysis
    end subroutine
end module