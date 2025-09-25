module transport_transient_m
    use diffusion_transient_m
    use transport_properties_heterog_m, only: tpt_props_heterog_c
    use transport_stab_params_m, only: stab_params_tpt_c
    implicit none
    save
    type, public, extends(diffusion_1D_transient_c) :: transport_1D_transient_c !> 1D transient transport subclass
        type(tpt_props_heterog_c) :: tpt_props_heterog      !> properties
        type(stab_params_tpt_c) :: stab_params_tpt          !> stability parameters
        logical :: Lagr_flag                                !> Lagrangian flag (if true, Lagrangian method is used, otherwise Eulerian method is used)
    contains
        procedure, public :: set_Lagr_flag
        procedure, public :: set_conc_r_flag=>set_conc_r_flag_tpt
        procedure, public :: compute_F_mat_PDE=>compute_F_mat_tpt
        procedure, public :: compute_trans_mat_PDE=>compute_trans_mat_tpt_transient
        procedure, public :: set_stab_params_tpt
        procedure, public :: initialise_PDE=>initialise_transport_1D_transient
        procedure, public :: initialise_transport_1D_transient_RT
        procedure, public :: write_PDE_1D=>write_transport_1D_transient
        procedure, public :: set_tpt_props_heterog_obj
        procedure, public :: mass_balance_error_ADE_trans_PMF_evap
        procedure, public :: mass_balance_error_ADE_trans_Dirichlet_evap
        procedure, public :: mass_balance_error_ADE_trans_PMF_recharge
        procedure, public :: mass_balance_error_ADE_trans_Dirichlet_recharge
        procedure, public :: mass_balance_error_ADE_trans_PMF_discharge
        procedure, public :: mass_balance_error_ADE_trans_Dirichlet_discharge
        procedure, public :: check_Delta_t
        procedure, public :: read_transport_data_WMA
        procedure, public :: write_transport_data_WMA
        procedure, public :: fund_sol_tpt_eqn_1D
        procedure, public :: compute_mixing_ratios_Delta_t_homog
        procedure, public :: update_mixing_ratios_Delta_t_homog
        procedure, public :: solve_tpt_EE_Delta_t_homog
    end type
    
    interface
        subroutine compute_F_mat_tpt(this)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c) :: this
        end subroutine
        
        subroutine compute_trans_mat_tpt_transient(this)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c) :: this
        end subroutine
        
        subroutine initialise_transport_1D_transient(this,root)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c) :: this
            character(len=*), intent(in) :: root !> root name for input files
        end subroutine
        
        subroutine initialise_transport_1D_transient_RT(this,root,mesh_type)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c) :: this
            !character(len=*), intent(in) :: path
            character(len=*), intent(in) :: root
            integer(kind=4), intent(in) :: mesh_type !> 1: 1D homogeneous, 2; 1D heterogeneous, 3: radial
            !character(len=*), intent(in) :: file_BCs
            !character(len=*), intent(in) :: file_spatial_discr
            !character(len=*), intent(in) :: file_time_discr
            !character(len=*), intent(in) :: file_tpt_props
        end subroutine
        
        subroutine write_transport_1D_transient(this,root,Time_out,output)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this !> transport object
            character(len=*), intent(in) :: root !> root name for output files
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(in) :: output(:,:)
        end subroutine
        
        function mass_balance_error_ADE_trans_Dirichlet_recharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_Dirichlet_discharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_Dirichlet_evap(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_PMF_recharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_PMF_discharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_PMF_evap(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        subroutine read_transport_data_WMA(this,root)!,mixing_ratios)!,f_vec)!,tpt_props,BCs,mesh,time_discr)
            import transport_1D_transient_c
            !import tpt_props_heterog_c
            !import BCs_t
            !import mesh_1D_Euler_homog_c
            !import time_discr_homog_c
            class(transport_1D_transient_c) :: this
            !character(len=*), intent(in) :: path
            !integer(kind=4), intent(in) :: unit
            character(len=*), intent(in) :: root
            !real(kind=8), intent(out), allocatable :: mixing_ratios(:,:)
            !real(kind=8), intent(out), allocatable :: f_vec(:)
            !type(tpt_props_heterog_c), intent(out) :: tpt_props
            !type(BCs_t), intent(out) :: BCs
            !!real(kind=8), intent(in) :: Delta_x
            !type(mesh_1D_Euler_homog_c), intent(out) :: mesh !> homogeneous Euler mesh 1D
            !type(time_discr_homog_c), intent(out) :: time_discr !> homogeneous time discretisation
        end subroutine
        
        subroutine read_discretisation_WMA(this,unit,root)
            import transport_1D_transient_c
            class(transport_1D_transient_c) :: this
            !character(len=*), intent(in) :: path
            integer(kind=4), intent(in) :: unit
            character(len=*), intent(in) :: root
        end subroutine
        
        subroutine write_transport_data_WMA(this,unit)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this                 !> 1D transient transport object                         
            integer(kind=4), intent(in) :: unit                                 !> unit of output file
        end subroutine
        
        subroutine compute_mixing_ratios_Delta_t_homog(this,A_mat_lumped)
            import transport_1D_transient_c
            import diag_matrix_c
            implicit none
            class(transport_1D_transient_c) :: this
            !real(kind=8), intent(in) :: theta
            type(diag_matrix_c), intent(out), optional :: A_mat_lumped
        end subroutine
        
        subroutine solve_tpt_EE_Delta_t_homog(this,Time_out,output)
            !> Solves 1D transient PDE with homogeneous time step using Lagr explicit method 
    
            !> this: transient PDE object
            !> Time_out: output time values
            !> output: concentration vs time output
    
            !> Results at all intermediate steps are written in binary mode in file conc_binary_EE.txt
    
            !use BCs_subroutines_m, only: Dirichlet_BCs_PDE, Neumann_homog_BCs, Robin_Neumann_homog_BCs
            import transport_1D_transient_c
    
            !> Variables
            class(transport_1D_transient_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        subroutine update_mixing_ratios_Delta_t_homog(this)
        import transport_1D_transient_c
        class(transport_1D_transient_c) :: this
        end subroutine
    end interface
    
    contains
      
        subroutine set_stab_params_tpt(this,stab_params_tpt)
            implicit none
            class(transport_1D_transient_c) :: this
            class(stab_params_tpt_c), intent(in) :: stab_params_tpt
            this%stab_params_tpt=stab_params_tpt
        end subroutine
        
        subroutine set_tpt_props_heterog_obj(this,tpt_props_heterog)
            implicit none
            class(transport_1D_transient_c) :: this
            class(tpt_props_heterog_c), intent(in) :: tpt_props_heterog
            this%tpt_props_heterog=tpt_props_heterog
        end subroutine
        
        subroutine set_conc_r_flag_tpt(this)
            implicit none
            class(transport_1D_transient_c) :: this
            
            integer(kind=4) :: i
            
            allocate(this%conc_obj%conc_r_flag(this%spatial_discr%Num_targets))
            this%conc_obj%conc_r_flag=0
            do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
                if (this%tpt_props_heterog%source_term(i)>0) then
                    this%conc_obj%conc_r_flag(i)=1
                end if
            end do
        end subroutine
        
        subroutine check_Delta_t(this)
            implicit none
            class(transport_1D_transient_c) :: this
            select type (time=>this%time_discr)
            type is (time_discr_homog_c)
                if (time%Delta_t>this%stab_params_tpt%Delta_t_crit) then
                    !print *, "Critical time step: ", this%stab_params_tpt%Delta_t_crit
                    !error stop "You must reduce time step to have stability"
                end if
            end select
        end subroutine

        function fund_sol_tpt_eqn_1D(this,c0,M,Delta_x,x,mu,t) result(conc)
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: c0
            real(kind=8), intent(in) :: M
            real(kind=8), intent(in) :: Delta_x
            real(kind=8), intent(in) :: x
            real(kind=8), intent(in) :: mu
            real(kind=8), intent(in) :: t
            !real(kind=8), intent(in) :: phi
            !real(kind=8), intent(in) :: D
            !real(kind=8), intent(in) :: q
            real(kind=8) :: conc
            
            real(kind=8), parameter :: eps=1d-12

                conc=exp(5d-1*(this%tpt_props_heterog%flux_cent(1)/this%tpt_props_heterog%dispersion(1))*&
                    (x-mu-t*this%tpt_props_heterog%flux_cent(1)/(2*this%tpt_props_heterog%porosity(1))))
                conc=conc*fund_sol_diff_eqn_1D(this,M,Delta_x,x,mu,t)+c0
        end function
        
        subroutine set_Lagr_flag(this,Lagr_flag)
        implicit none
        class(transport_1D_transient_c) :: this
        logical, intent(in) :: Lagr_flag
        this%Lagr_flag=Lagr_flag
        end subroutine
end module