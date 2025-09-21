!> Boundary conditions type:
!>   contains information of boundary conditions
module BCs_m
    use time_fct_m !> time function module
    implicit none
    save
    type, public :: BCs_t
        integer(kind=4) :: BCs_label(2)         !> First element: inflow
                                                !> second element: outflow
                                                !> 1: Dirichlet
                                                !> 2: Neumann homogeneous
                                                !> 3: Robin
        logical :: evap                         !> evaporation flag
        real(kind=8) :: conc_inf                !> concentration at inflow
        real(kind=8) :: conc_out                !> concentration at outflow
        real(kind=8) :: flux_inf                !> flux at inflow
        real(kind=8) :: flux_out                !> flux at outflow
        real(kind=8) :: head_inf                !> head at inflow
        real(kind=8) :: head_inf_D              !> dimensionless head at inflow
        real(kind=8) :: head_out                !> head at outflow
        real(kind=8) :: head_out_D              !> dimensionless head at outflow
        real(kind=8) :: caudal_inf              !> caudal at inflow
        type(time_fct_real_c) :: flow_inf       !> caudal (Q) time function
    contains
        procedure, public :: set_BCs_label
        procedure, public :: set_evap
        procedure, public :: set_flow_inf
        procedure, public :: read_BCs
        procedure, public :: read_Dirichlet_BCs_conc
        procedure, public :: read_Dirichlet_BCs_head
        procedure, public :: read_Robin_BC_inflow
        procedure, public :: read_caudal_inf
        procedure, public :: set_conc_boundary
        procedure, public :: set_cst_flux_boundary
        procedure, public :: compute_dimless_BCs
    end type
    
    contains
        subroutine set_BCs_label(this,BCs)
            implicit none
            class(BCs_t) :: this
            integer(kind=4), intent(in) :: BCs(2)
            if (BCs(1)>3 .or. BCs(2)>3) error stop "BCs not implemented yet"
            this%BCs_label=BCs
        end subroutine
        
        subroutine set_evap(this,evap)
            implicit none
            class(BCs_t) :: this
            logical, intent(in) :: evap
            this%evap=evap
        end subroutine 
        
        subroutine set_conc_boundary(this,conc_inf,conc_out)
            implicit none
            class(BCs_t) :: this
            real(kind=8), intent(in) :: conc_inf,conc_out
            this%conc_inf=conc_inf
            this%conc_out=conc_out
        end subroutine
        
        subroutine set_cst_flux_boundary(this,flux)
            implicit none
            class(BCs_t) :: this
            real(kind=8), intent(in) :: flux
            this%flux_inf=flux
            this%flux_out=this%flux_inf
        end subroutine
        
        subroutine read_BCs(this,filename)
            implicit none
            class(BCs_t) :: this
            character(len=*), intent(in) :: filename
            open(unit=1,file=filename,status='old',action='read')
            read(1,*) this%BCs_label
            read(1,*) this%evap
            close(1)
        end subroutine
        
        subroutine read_Dirichlet_BCs_conc(this,filename)
            implicit none
            class(BCs_t) :: this
            character(len=*), intent(in) :: filename
            open(unit=1,file=filename,status='old',action='read')
            read(1,*) this%conc_inf
            read(1,*) this%conc_out
            close(1)
        end subroutine
        
        subroutine read_Robin_BC_inflow(this,filename)
            implicit none
            class(BCs_t) :: this
            character(len=*), intent(in) :: filename
            open(unit=1,file=filename,status='old',action='read')
            !read(1,*) this%flux_inf
            read(1,*) this%caudal_inf !> caudal at inflow
            read(1,*) this%conc_inf !> concentration at inflow (not used in WMA)
            close(1)
        end subroutine
        
        subroutine read_caudal_inf(this,filename)
            implicit none
            class(BCs_t) :: this
            character(len=*), intent(in) :: filename
            open(unit=1,file=filename,status='old',action='read')
            read(1,*) this%caudal_inf
            close(1)
            !this%flux_inf=this%caudal_inf
        end subroutine
        
        subroutine set_flow_inf(this,flow_inf)
        implicit none
        class(BCs_t) :: this
        type(time_fct_real_c), intent(in) :: flow_inf
        this%flow_inf=flow_inf
        end subroutine
        
        subroutine read_Dirichlet_BCs_head(this,filename)
            implicit none
            class(BCs_t) :: this
            character(len=*), intent(in) :: filename
            open(unit=1,file=filename,status='old',action='read')
            read(1,*) this%head_inf
            read(1,*) this%head_out
            close(1)
        end subroutine
        
        subroutine compute_dimless_BCs(this,h_c)
            implicit none
            class(BCs_t) :: this
            real(kind=8), intent(in) :: h_c !> characteristic head
            this%head_inf_D=this%head_inf/h_c
            this%head_out_D=this%head_out/h_c
        end subroutine
end module