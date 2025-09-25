module flow_props_heterog_m
    use properties_m
    implicit none
    save
    type, public, extends(props_c) :: flow_props_heterog_c !> heterogeneous flow properties subclass
        !> physical properties
        real(kind=8), allocatable :: spec_stor(:) !> specific storage (S_s)
        real(kind=8), allocatable :: hydr_cond(:) !> hydraulic conductivity (K)
    contains
        procedure, public :: set_props_flow_heterog
        procedure, public :: read_props=>read_props_flow_heterog
        procedure, public :: are_props_homog=>are_flow_props_homog
    end type
    
    type, public, extends(flow_props_heterog_c) :: flow_props_heterog_conf_c !> heterogeneous 1D diffusion properties subclass
        !> physical properties
        real(kind=8), allocatable :: storativity(:)   !> S
        real(kind=8), allocatable :: transmissivity(:) !> T
    contains
        procedure, public :: read_props=>read_props_flow_heterog_conf
    end type
    
    contains
        subroutine set_props_flow_heterog(this,porosity,dispersion)
            implicit none
            class(flow_props_heterog_c) :: this
            real(kind=8), intent(in) :: porosity(:),dispersion(:)
            !if (size(porosity)/=size(dispersion)) error stop "Dimensions of porosity and dispersion must be the same"
            !this%porosity=porosity
            !this%dispersion=dispersion
        end subroutine
        
        subroutine read_props_flow_heterog(this,root,spatial_discr)
            implicit none
            class(flow_props_heterog_c) :: this
            character(len=*), intent(in) :: root
            class(spatial_discr_c), intent(in), optional :: spatial_discr
            
            real(kind=8), parameter :: epsilon=1d-12
            real(kind=8) :: S_s,K,w
            logical :: flag
            
            call this%set_homog_flag(.true.) !> by default
            call this%set_stat_flag(.true.) !> by default
            
            open(unit=1,file=root//'_flow_props.dat',status='old',action='read')
            read(1,*) flag
            if (flag.eqv..true.) then
                backspace(1)
                read(1,*) flag, w
                allocate(this%source_term(spatial_discr%Num_targets))
                this%source_term=w
                this%source_term_order=0 !> constant source term
                !if (abs(r)<epsilon) then
                !    this%source_term_order=0
                !else
                !    this%source_term_order=1
                !end if
            !else if (allocated(this%flux)) then
            !    continue
            else
                !> Chapuza
                !allocate(this%source_term(spatial_discr%Num_targets))
                !read(1,*) this%source_term(1)
                !this%source_term(2:)=this%source_term(1)
                !this%source_term_order=1
                !allocate(this%source_term(spatial_discr%Num_targets))
                this%homog_flag=.false.
                read(1,*) this%source_term
                if (size(this%source_term)/=spatial_discr%Num_targets) error stop "Dimension error in source term"
            end if
            read(1,*) flag
            if (flag .eqv. .true.) then
                backspace(1)
                read(1,*) flag, S_s
                allocate(this%spec_stor(spatial_discr%Num_targets-spatial_discr%targets_flag))
                this%spec_stor=S_s
            else
                this%homog_flag=.false.
            end if
            read(1,*) flag
            if (flag.eqv..true.) then
                backspace(1)
                read(1,*) flag, K
                allocate(this%hydr_cond(spatial_discr%Num_targets))
                this%hydr_cond=K
            else
                this%homog_flag=.false.
            end if
            !read(1,*) flag !> flux flag
            !if ((flag .eqv. .true.) .and. (allocated(this%flux).eqv..true.)) then
            !    error stop "Flux already allocated"
            !else if (flag.eqv..true. .and. sum_squares(this%source_term)>epsilon) then
            !    error stop "Flux cannot be constant"
            !else if (flag.eqv..true.) then
            !    this%cst_flux_flag=.true. !> constant flux
            !    backspace(1) !> chapuza
            !    read(1,*) flag, q !> we read the flux
            !    if ((spatial_discr%scheme == 2) .and. (spatial_discr%targets_flag == 0)) then
            !        n_flux=spatial_discr%Num_targets+1
            !    else
            !        n_flux=spatial_discr%Num_targets
            !    end if
            !    allocate(this%flux(n_flux))
            !    this%flux=q
            !else
            !    this%homog_flag=.false.
            !end if
            close(1)
        end subroutine
        
        subroutine read_props_flow_heterog_conf(this,root,spatial_discr)
            implicit none
            class(flow_props_heterog_conf_c) :: this
            character(len=*), intent(in) :: root
            class(spatial_discr_c), intent(in), optional :: spatial_discr
            
            real(kind=8), parameter :: epsilon=1d-12
            real(kind=8) :: S,T,r
            logical :: flag
            
            call this%set_homog_flag(.true.) !> by default
            call this%set_stat_flag(.true.) !> by default
            
            open(unit=1,file=root//'_flow_props.dat',status='old',action='read')
            read(1,*) flag
            if (flag.eqv..true.) then
                backspace(1)
                read(1,*) flag, r
                allocate(this%source_term(spatial_discr%Num_targets))
                this%source_term=r
                this%source_term_order=0 !> constant source term
                !if (abs(r)<epsilon) then
                !    this%source_term_order=0
                !else
                !    this%source_term_order=1
                !end if
            !else if (allocated(this%flux)) then
            !    continue
            else
                !> Chapuza
                !allocate(this%source_term(spatial_discr%Num_targets))
                !read(1,*) this%source_term(1)
                !this%source_term(2:)=this%source_term(1)
                !this%source_term_order=1
                !allocate(this%source_term(spatial_discr%Num_targets))
                this%homog_flag=.false.
                read(1,*) this%source_term
                if (size(this%source_term)/=spatial_discr%Num_targets) error stop "Dimension error in source term"
            end if
            read(1,*) flag
            if (flag.eqv..true.) then
                backspace(1)
                read(1,*) flag, S
                allocate(this%storativity(spatial_discr%Num_targets-spatial_discr%targets_flag))
                this%storativity=S
            else
                this%homog_flag=.false.
            end if
            read(1,*) flag
            if (flag.eqv..true.) then
                backspace(1)
                read(1,*) flag, T
                allocate(this%transmissivity(spatial_discr%Num_targets+1-spatial_discr%targets_flag))
                this%transmissivity=T
            else 
                this%homog_flag=.false.
            end if
            !read(1,*) flag !> flux flag
            !if ((flag .eqv. .true.) .and. (allocated(this%flux).eqv..true.)) then
            !    error stop "Flux already allocated"
            !else if (flag.eqv..true. .and. sum_squares(this%source_term)>epsilon) then
            !    error stop "Flux cannot be constant"
            !else if (flag.eqv..true.) then
            !    this%cst_flux_flag=.true. !> constant flux
            !    backspace(1) !> chapuza
            !    read(1,*) flag, q !> we read the flux
            !    if ((spatial_discr%scheme == 2) .and. (spatial_discr%targets_flag == 0)) then
            !        n_flux=spatial_discr%Num_targets+1
            !    else
            !        n_flux=spatial_discr%Num_targets
            !    end if
            !    allocate(this%flux(n_flux))
            !    this%flux=q
            !else
            !    this%homog_flag=.false.
            !end if
            close(1)
        end subroutine
        
      subroutine are_flow_props_homog(this)
            implicit none
            class(flow_props_heterog_c) :: this
            
            integer(kind=4) :: i
            real(kind=8), parameter :: eps=1d-12
            
            this%homog_flag=.true.
            !do i=2,size(this%porosity) !> we assume porosity & dispersion have the same dimension
            !    if (abs(this%dispersion(1)-this%dispersion(i))>eps .or. abs(this%porosity(1)-this%porosity(i))>eps) then
            !        this%homog_flag=.false.
            !        exit
            !    end if
            !end do
        end subroutine
        
        
end module