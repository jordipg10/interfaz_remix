!> 1D steady-state transport equation:
!> 0 = -q*c' + D*c'' + r*(c_r-c)
module transport_m
    use diffusion_m, only: diffusion_1D_c, PDE_1D_c
    use transport_properties_heterog_m, only: tpt_props_heterog_c
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c, spatial_discr_c
    implicit none
    save
    type, public, extends(diffusion_1D_c) :: transport_1D_c
        type(tpt_props_heterog_c) :: tpt_props_heterog
    contains
    !> set
        procedure, public :: set_tpt_props_heterog_obj
        procedure, public :: set_conc_r_flag=>set_conc_r_flag_tpt
    !> Computations
        procedure, public :: compute_trans_mat_PDE=>compute_trans_mat_tpt
        procedure, public :: mass_balance_error_ADE_stat_Dirichlet_discharge
    !> Abstract
        procedure, public :: initialise_PDE=>initialise_transport_1D
        procedure, public :: write_PDE_1D=>write_transport_1D
    end type
!**************************************************************************************************
    interface
        subroutine compute_trans_mat_tpt(this)
            import transport_1D_c
            implicit none
            class(transport_1D_c) :: this
        end subroutine
        
        subroutine initialise_transport_1D(this,root)
            import transport_1D_c
            implicit none
            class(transport_1D_c) :: this
            character(len=*), intent(in) :: root !> root name for input files
        end subroutine
        
        subroutine write_transport_1D(this,root,Time_out,output)
            import transport_1D_c
            implicit none
            class(transport_1D_c), intent(in) :: this
            character(len=*), intent(in) :: root !> root name for output files
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(in) :: output(:,:)
        end subroutine
        
        function mass_balance_error_ADE_stat_Dirichlet_discharge(this) result(mass_bal_err)
            import transport_1D_c
            implicit none
            class(transport_1D_c), intent(in) :: this
            real(kind=8) :: mass_bal_err
        end function
    
    end interface
!*****************************************************************************************************************************
    contains
        
        subroutine set_tpt_props_heterog_obj(this,tpt_props_heterog)
            implicit none
            class(transport_1D_c) :: this
            class(tpt_props_heterog_c), intent(in)  :: tpt_props_heterog
            this%tpt_props_heterog=tpt_props_heterog
        end subroutine
        
        subroutine set_conc_r_flag_tpt(this)
            implicit none
            class(transport_1D_c) :: this
            
            integer(kind=4) :: i
            
            allocate(this%conc_r_flag(this%spatial_discr%Num_targets))
            this%conc_r_flag=0
            do i=1,this%spatial_discr%Num_targets
                if (this%tpt_props_heterog%source_term(i)>0) then
                    this%conc_r_flag(i)=1
                end if
            end do
        end subroutine 

        ! function anal_sol_tpt_1D_stat_flujo_lin(this,x) result(conc)
        !     !> Analytical solution 1D stationary transport equation with:
        !     !>   D=cst
        !     !>   q(x)=-x
        !     !>   c(0)=1
        !     !>   c(L)=0
        !         class(transport_1D_c), intent(in) :: this
        !         real(kind=8), intent(in) :: x
        !         real(kind=8) :: conc
                
        !         real(kind=8) :: D
        !         real(kind=8), parameter :: pi=4d0*atan(1d0)
                
        !         D=this%tpt_props_heterog%dispersion(1)
    
        !         if (this%dimensionless.eqv..true.) then
        !             conc=-erf(x/sqrt(2d0))/erf(1d0/sqrt(2d0)) + 1
        !         else
        !             conc=-erf(x/sqrt(2d0*D))/erf(this%spatial_discr%measure/sqrt(2d0*D)) + 1d0
        !         end if
        !     end function
            
       
        !     function der_anal_sol_tpt_1D_stat_flujo_lin(this,x) result(der_conc)
        !     !> Derivative of analytical solution 1D stationary transport equation with:
        !     !>   D=cst
        !     !>   q(x)=-x
        !     !>   c(0)=1
        !     !>   c(L)=0
        !         class(transport_1D_c), intent(in) :: this
        !         real(kind=8), intent(in) :: x
        !         real(kind=8) :: der_conc
                
        !         real(kind=8) :: D,L 
        !         real(kind=8), parameter :: pi=4d0*atan(1d0)
                
        !         L=this%spatial_discr%measure !> 1D
        !         D=this%tpt_props_heterog%dispersion(1)
        !         der_conc=(-sqrt(2d0)/(erf(this%spatial_discr%measure/sqrt(2d0*D))*sqrt(pi*D)))*exp(-(x**2)/(2d0*D)) !> derivative analytical solution
        !     end function
    
            
        
        !     function der_anal_sol_tpt_1D_stat_flujo_cuad(this,x) result(der_conc)
        !     !> Derivative of analytical solution 1D stationary transport equation with:
        !     !>   D=cst
        !     !>   q(x)=-x^2
        !     !>   c(0)=1
        !     !>   c(L)=0
        !         class(transport_1D_c), intent(in) :: this
        !         real(kind=8), intent(in) :: x
        !         real(kind=8) :: der_conc
                
        !         real(kind=8) :: D,L
        !         real(kind=8), parameter :: pi=4d0*atan(1d0),incompl_gamma_term=0.327336564991358
                
        !         D=this%tpt_props_heterog%dispersion(1)
        !         L=this%spatial_discr%measure !> 1D
    
        !         der_conc=-exp(-(x**3)/(3d0*D))*(3d0**(2d0/3d0))/(gamma(1d0/3d0)-incompl_gamma_term)
        !     end function
        
end module 