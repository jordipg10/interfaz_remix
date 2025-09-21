!> 1D steady-state diffusion equation:
!> 0 = D*c'' + r*(c_r-c)
module diffusion_m
    use PDE_m, only: PDE_1D_c
    use diff_props_heterog_m, only: diff_props_heterog_c, props_c
    implicit none
    save
    type, public, extends(PDE_1D_c) :: diffusion_1D_c !> 1D diffusion equation class
        real(kind=8), allocatable :: conc(:) !> concentrations (c)
        real(kind=8), allocatable :: conc_ext(:) !> (c_ext) external concentrations (include recharge & boundary)
        real(kind=8), allocatable :: conc_r(:) !> (c_r) recharge concentrations (c_r>0)
        real(kind=8), allocatable :: conc_bd(:) !> (c_bd) boundary concentrations (c_bd>0)
        integer(kind=4), allocatable :: conc_r_flag(:)      !> 1 if r>0
                                                            !> 0 otherwise
        type(diff_props_heterog_c) :: diff_props_heterog
    contains
    !> set
        procedure, public :: set_conc_ext
        procedure, public :: set_diff_props_heterog
        procedure, public :: set_conc_r_flag=>set_conc_r_flag_diff
    !> Computations
        procedure, public :: compute_trans_mat_PDE=>compute_trans_mat_diff
        procedure, public :: compute_source_term_PDE=>compute_source_term_diff
        procedure, public :: compute_rech_mat_PDE=>compute_rech_mat_diff
        procedure, public :: initialise_PDE=>initialise_diffusion_1D
        procedure, public :: write_PDE_1D=>write_diffusion_1D
        procedure, public :: allocate_conc
        procedure, public :: solve_PDE_1D=>solve_diff_1D
    end type
!*****************************************************************************************************************************
    interface
        
        subroutine initialise_diffusion_1D(this,root)
            import diffusion_1D_c
            class(diffusion_1D_c) :: this
            character(len=*), intent(in) :: root !> root name for input files
        end subroutine
        
        subroutine compute_trans_mat_diff(this)
            import diffusion_1D_c
            implicit none
            class(diffusion_1D_c) :: this
        end subroutine
        
       subroutine write_diffusion_1D(this,root,Time_out,output)
            import diffusion_1D_c
            import props_c
            class(diffusion_1D_c), intent(in) :: this
            character(len=*), intent(in) :: root !> root name for output files
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(in) :: output(:,:)
        end subroutine
      
        subroutine solve_diff_1D(this,Time_out,output)
        import diffusion_1D_c
        class(diffusion_1D_c) :: this
        real(kind=8), intent(in) :: Time_out(:)
        real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        
        
    end interface
!*****************************************************************************************************************************
    contains
        
        subroutine set_conc_ext(this,conc_ext)
            class(diffusion_1D_c) :: this
            real(kind=8), intent(in) :: conc_ext(:)
            if (size(conc_ext)/=this%spatial_discr%Num_targets) error stop "Dimension error in external concentration"
            this%conc_ext=conc_ext
        end subroutine 
        
        subroutine set_diff_props_heterog(this,diff_props_heterog)
            implicit none
            class(diffusion_1D_c) :: this
            class(diff_props_heterog_c), intent(in) :: diff_props_heterog
            this%diff_props_heterog=diff_props_heterog
        end subroutine
        
        subroutine set_conc_r_flag_diff(this)
            implicit none
            class(diffusion_1D_c) :: this
            
            integer(kind=4) :: i
            
            allocate(this%conc_r_flag(this%spatial_discr%Num_targets))
            this%conc_r_flag=0
            do i=1,this%spatial_discr%Num_targets
                if (this%diff_props_heterog%source_term(i)>0) then
                    this%conc_r_flag(i)=1
                end if
            end do
        end subroutine 
        
        subroutine compute_source_term_diff(this)
            !> $g=E*c_ext$
            !use PDE_m, only: PDE_1D_c
            !use transport_m, only: transport_1D_c, diffusion_1D_c
            !use transport_transient_m, only: transport_1D_transient_c, diffusion_1D_transient_c
            implicit none
            class(diffusion_1D_c) :: this
    
            !allocate(this%source_term_PDE(this%spatial_discr%Num_targets))
            !this%source_term_PDE=0d0 !> $g=0$ chapuza
            !select type (this)
            !class is (diffusion_1D_c)
            !    this%source_term_PDE=this%rech_mat%diag*this%conc_ext
            !class is (diffusion_1D_transient_c)
            !    this%source_term_PDE=this%rech_mat%diag*this%conc_ext
            !end select
            this%source_term_PDE=this%rech_mat%diag*this%conc_ext
            this%source_term_PDE(1)=this%bd_mat(1)*this%BCs%conc_inf
            this%source_term_PDE(this%spatial_discr%Num_targets)=this%bd_mat(2)*this%BCs%conc_out
            ! select type (this)
            ! type is (transport_1D_c)
            !     this%source_term_PDE=this%conc_r_flag*this%tpt_props_heterog%source_term*this%conc_ext
            ! type is (transport_1D_transient_c)
            !     this%source_term_PDE=this%conc_r_flag*this%tpt_props_heterog%source_term*this%conc_ext
            ! type is (diffusion_1D_c)
            !     this%source_term_PDE=this%diff_props_heterog%source_term*this%conc_ext
            ! type is (diffusion_1D_transient_c)
            !     this%source_term_PDE=this%diff_props_heterog%source_term*this%conc_ext
            ! end select

        end subroutine
        
        subroutine allocate_conc(this)
            implicit none
            class(diffusion_1D_c) :: this
            allocate(this%conc(this%spatial_discr%Num_targets), this%conc_ext(this%spatial_discr%Num_targets), &
                this%conc_r(this%spatial_discr%Num_targets), this%conc_bd(this%spatial_discr%Num_targets), &
                this%conc_r_flag(this%spatial_discr%Num_targets))
        end subroutine
        
        subroutine compute_rech_mat_diff(this)
        implicit none
        class(diffusion_1D_c) :: this
        !> $R$ matrix for recharge
        !> $R=diag(r)$
        this%rech_mat%diag=this%diff_props_heterog%source_term
        end subroutine
end module 