module diff_props_heterog_m
    use properties_m
    implicit none
    save
    type, public, extends(props_c) :: diff_props_heterog_c !> heterogeneous 1D diffusion properties subclass
        !> physical properties
        real(kind=8), allocatable :: porosity(:)   !> phi
        real(kind=8), allocatable :: dispersion(:) !> D
        ! real(kind=8) :: long_disp !> longitudinal dispersivity (alpha_L)
    contains
        procedure, public :: set_props_diff_heterog
        procedure, public :: read_props=>read_props_diff_heterog
        procedure, public :: are_props_homog=>are_diff_props_homog
    end type
    
    contains
        subroutine set_props_diff_heterog(this,porosity,dispersion)
            implicit none
            class(diff_props_heterog_c) :: this
            real(kind=8), intent(in) :: porosity(:),dispersion(:)
            if (size(porosity)/=size(dispersion)) error stop "Dimensions of porosity and dispersion must be the same"
            this%porosity=porosity
            this%dispersion=dispersion
        end subroutine
        
        subroutine read_props_diff_heterog(this,root,spatial_discr)
            implicit none
            class(diff_props_heterog_c) :: this
            character(len=*), intent(in) :: root
            class(spatial_discr_c), intent(in), optional :: spatial_discr
            
            real(kind=8), parameter :: epsilon=1d-12
            real(kind=8) :: phi,D,r
            integer(kind=4) :: flag
            
            open(unit=1,file=root//'_diff_props.dat',status='old',action='read')
            read(1,"(/,F10.2)") flag
            if (flag == 1) then
                backspace(1)
                read(1,*) flag, r
                allocate(this%source_term(spatial_discr%Num_targets-spatial_discr%targets_flag))
                this%source_term=r
                this%source_term_order=0
            else
                read(1,*) this%source_term
                if (size(this%source_term)/=spatial_discr%Num_targets-spatial_discr%targets_flag) error stop "Dimension error in &
                    source term"
            end if
            read(1,*) flag
            if (flag == 1) then
                backspace(1)
                read(1,*) flag, phi
                allocate(this%porosity(spatial_discr%Num_targets-spatial_discr%targets_flag))
                this%porosity=phi
            end if
            read(1,*) flag
            if (flag == 1) then
                backspace(1)
                read(1,*) flag, D
                allocate(this%dispersion(spatial_discr%Num_targets-spatial_discr%targets_flag))
                this%dispersion=D
            end if
            close(1)
        end subroutine
        
      subroutine are_diff_props_homog(this)
            implicit none
            class(diff_props_heterog_c) :: this
            
            integer(kind=4) :: i
            real(kind=8), parameter :: eps=1d-12
            
            this%homog_flag=.true.
            do i=2,size(this%porosity) !> we assume porosity & dispersion have the same dimension
                if (abs(this%dispersion(1)-this%dispersion(i))>eps .or. abs(this%porosity(1)-this%porosity(i))>eps) then
                    this%homog_flag=.false.
                    exit
                end if
            end do
        end subroutine
        
        
end module