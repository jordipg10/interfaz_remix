module transport_properties_heterog_m
    use diff_props_heterog_m, only: diff_props_heterog_c, are_diff_props_homog
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c, spatial_discr_c
    !use diffusion_m, only: diffusion_1D_c, PDE_1D_c
    !use transport_m, only: transport_1D_c
    !use transport_transient_m, only: transport_1D_transient_c, diffusion_1D_transient_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    use polynomials_m, only: real_poly_1D, der_real_poly_1D, sec_der_real_poly_1D
    use vectors_m, only: sum_squares
    use time_fct_m, only: time_fct_real_c
    implicit none
    save
    type, public, extends(diff_props_heterog_c) :: tpt_props_heterog_c !> heterogeneous 1D transport properties subclass
        real(kind=8), allocatable :: flux_int(:) !> flux at the interfaces
        real(kind=8), allocatable :: flux_cent(:) !> flux at the centre of the cells/elements
        real(kind=8), allocatable :: dispersivity(:) !> alpha
        !real(kind=8), allocatable :: dispersion(:) !> D
        logical :: cst_flux_flag !> TRUE if flux is constant, FALSE otherwise
        !integer(kind=4) :: num_flux_chg !> number of time steps where flux changes
        !integer(kind=4), allocatable :: flux_change(:) !> array of time steps where flux changes
        type(time_fct_real_c) :: flux_trans !> transient flux (we assume homogeneous)
    contains
        procedure, public :: set_cst_flux_flag !> set constant flux flag
        procedure, public :: update_flux_int_trans
        procedure, public :: read_props=>read_tpt_props_heterog
        procedure, public :: compute_dispersion
        procedure, public :: compute_flux_lin
        procedure, public :: compute_flux_nonlin
        procedure, public :: compute_flux_rad
        procedure, public :: compute_source_term
        procedure, public :: are_props_homog=>are_tpt_props_homog
        procedure, public :: allocate_flux
    end type
    
    contains
        subroutine read_source_term_tpt(this,filename,mesh)
            implicit none
            class(tpt_props_heterog_c) :: this
            character(len=*), intent(in) :: filename
            class(spatial_discr_c), intent(in) :: mesh
            
            real(kind=8) :: r
            logical :: cst_source_term
            
            open(unit=1,file=filename,status='old',action='read')
            read(1,*) cst_source_term
            if (cst_source_term.eqv..true.) then
                backspace(1)
                read(1,*) cst_source_term, r
                allocate(this%source_term(mesh%Num_targets-mesh%targets_flag))
                this%source_term=r
                this%source_term_order=0
            else if (allocated(this%flux_int) .and. allocated(this%flux_cent)) then
                continue
            else
                read(1,*) this%source_term
                if (size(this%source_term)/=mesh%Num_targets-mesh%targets_flag) error stop "Dimension error in source term"
            end if
            close(1)
        end subroutine
        
        subroutine set_tpt_props_heterog(this,porosity,dispersion,flux_int,flux_cent)
            implicit none
            class(tpt_props_heterog_c) :: this
            real(kind=8), intent(in) :: porosity(:)
            real(kind=8), intent(in) :: dispersion(:)
            real(kind=8), intent(in) :: flux_int(:)
            real(kind=8), intent(in) :: flux_cent(:)
            if (size(flux_cent)/=size(porosity)) then
                error stop "Dimensions of porosity and flux at centres must be the same"
            else if (size(flux_int)/=size(flux_cent)+1) then
                error stop "Dimensions of fluxes are wrong"
            else
                this%porosity=porosity
                this%dispersion=dispersion
                this%flux_int=flux_int
                this%flux_cent=flux_cent
            end if
        end subroutine
        
        subroutine read_tpt_props_heterog(this,root,spatial_discr)
            implicit none
            class(tpt_props_heterog_c) :: this
            character(len=*), intent(in) :: root
            class(spatial_discr_c), intent(in), optional :: spatial_discr
            
            integer(kind=4) :: n_flux,Num_cells
            real(kind=8), parameter :: epsilon=1d-12
            real(kind=8) :: phi,D,q,r,alpha
            logical :: flag,trans_flag
            
            call this%set_homog_flag(.true.) !> by default
            call this%set_stat_flag(.true.) !> by default
            
            open(unit=1,file=root//'_tpt_props.dat',status='old',action='read')
            read(1,*) flag !> source term flag
            allocate(this%source_term(spatial_discr%Num_targets))
            if (flag.eqv..true.) then
                backspace(1)
                read(1,*) flag, r
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
                call this%set_homog_flag(.false.)
                !read(1,*) this%source_term
                !if (size(this%source_term)/=spatial_discr%Num_targets) error stop "Dimension error in source term"
            end if
            read(1,*) flag !> porosity flag
            if (flag.eqv..true.) then
                backspace(1)
                read(1,*) flag, phi
                allocate(this%porosity(spatial_discr%Num_targets))
                this%porosity=phi
            else
                allocate(this%porosity(spatial_discr%Num_targets))
                !this%porosity=1d0 !> default porosity
            end if
            read(1,*) flag !> constant flux flag
            if ((flag .eqv. .true.) .and. (allocated(this%flux_int).eqv..true.)) then
                error stop "Flux already allocated"
            ! else if (flag.eqv..true. .and. sum_squares(this%source_term)>epsilon) then
            !     error stop "Flux cannot be constant"
            else! if (flag.eqv..true.) then
                call this%set_cst_flux_flag(flag) !> constant flux flag
                !if ((spatial_discr%scheme == 2) .and. (spatial_discr%targets_flag == 0)) then
                !    n_flux=spatial_discr%Num_targets+1
                !else
                !    n_flux=spatial_discr%Num_targets
                !end if
                Num_cells=spatial_discr%Num_targets-spatial_discr%targets_flag
                allocate(this%flux_cent(Num_cells),this%flux_int(Num_cells+1))
                if (flag.eqv..false.) then !> non-constant flux
                    call this%set_homog_flag(flag) !> homogeneous flag
                    backspace(1) !> chapuza
                    read(1,*) flag, trans_flag !> we read the transient flag
                else !> constant flux
                    backspace(1) !> chapuza
                    read(1,*) flag, q, trans_flag !> we read the flux and the transient flag
                    this%flux_int=q
                    this%flux_cent=q
                end if
                if (trans_flag.eqv..true.) then
                    call this%set_stat_flag(.false.)
                    open(unit=60,file=root//'_flow_inf.dat',status='old',action='read')
                    call this%flux_trans%read_time_series(60) !> we assume homogeneous flux
                    close(60)
                    !this%flux_trans%set_time(spatial_discr%Delta_t)
                end if
            ! else
            !     call this%set_homog_flag(.false.)
            end if
            read(1,*) flag !> dispersivity flag
            if (flag.eqv..true.) then
                backspace(1)
                read(1,*) flag, alpha
                !> dispersivity dimension depends on scheme and mesh type
                if (spatial_discr%scheme == 2 .and. spatial_discr%targets_flag == 0) then !> cell centred, proposed by Petcham� & Carrera (2024)
                    allocate(this%dispersivity(spatial_discr%Num_targets+1))
                !else if (spatial_discr%scheme == 1 .and. spatial_discr%targets_flag == 0) then !> cell centred, traditional
                !    allocate(this%dispersivity(spatial_discr%Num_targets))
                else
                    allocate(this%dispersivity(spatial_discr%Num_targets))
                end if
                this%dispersivity=alpha
            end if
            close(1)
        end subroutine
        
        subroutine compute_flux_lin(this,q_inf,spatial_discr_obj)
            implicit none
            class(tpt_props_heterog_c) :: this
            real(kind=8), intent(in) :: q_inf !> flux entering domain
            class(spatial_discr_c), intent(in) :: spatial_discr_obj
            !real(kind=8), intent(out) :: q_out !> flux leaving domain
            
            integer(kind=4) :: Num_cells,i
            
            Num_cells=spatial_discr_obj%Num_targets-spatial_discr_obj%targets_flag
            
            !if ((spatial_discr_obj%scheme == 2) .and. (spatial_discr_obj%targets_flag == 0)) then !> IFDS
                call this%allocate_flux(Num_cells+1)
                this%flux_int(1)=q_inf
                select type (spatial_discr_obj)
                type is (mesh_1D_Euler_homog_c)
                    do i=1,Num_cells
                        this%flux_int(i+1)=this%flux_int(i)+this%source_term(i)*spatial_discr_obj%Delta_x
                        this%flux_cent(i)=5d-1*(this%flux_int(i)+this%flux_int(i+1)) !> average at cell centre
                    end do
                type is (mesh_1D_Euler_heterog_c)
                    do i=1,Num_cells
                        this%flux_int(i+1)=this%flux_int(i)+this%source_term(i)*spatial_discr_obj%Delta_x(i)
                        this%flux_cent(i)=5d-1*(this%flux_int(i)+this%flux_int(i+1)) !> average at cell centre
                    end do
                end select
            !    q_out=this%flux(Num_cells+1)
            !else if (spatial_discr_obj%scheme == 1 .and. spatial_discr_obj%targets_flag == 0) then
            !    allocate(this%flux(Num_cells))
            !    select type (spatial_discr_obj)
            !    type is (mesh_1D_Euler_homog_c)
            !        this%flux(1)=q_inf+this%source_term(1)*spatial_discr_obj%Delta_x/2d0
            !        do i=2,Num_cells
            !            this%flux(i)=this%flux(i-1)+this%source_term(i-1)*spatial_discr_obj%Delta_x
            !        end do
            !    type is (mesh_1D_Euler_heterog_c)
            !        this%flux(1)=q_inf+this%source_term(1)*spatial_discr_obj%Delta_x(1)/2d0
            !        do i=2,Num_cells
            !            this%flux(i)=this%flux(i-1)+this%source_term(i-1)*(spatial_discr_obj%Delta_x(i-1)+&
            !            spatial_discr_obj%Delta_x(i))/2d0
            !        end do
            !    end select
            !    q_out=this%flux(Num_cells)+this%source_term(Num_cells)*spatial_discr_obj%get_Cell_size(Num_cells)/2d0
            !end if
        end subroutine
        
        subroutine compute_flux_nonlin(this,flux_coeffs,spatial_discr_obj) !> we assume domain starts at x=0
            implicit none
            class(tpt_props_heterog_c) :: this
            real(kind=8), intent(in) :: flux_coeffs(:) !> coefficients of flux polynomial in decreasing order
            class(spatial_discr_c), intent(in) :: spatial_discr_obj
            !real(kind=8), intent(out) :: q_out !> flux leaving domain
            
            integer(kind=4) :: i,deg,Num_cells
            real(kind=8) :: x,dqn_dx,d2qn_dx2,Delta_x_n,x0
            
            deg=size(flux_coeffs)-1
            Num_cells=spatial_discr_obj%Num_targets-spatial_discr_obj%targets_flag
            
            !if (spatial_discr_obj%scheme == 1) then
                allocate(this%flux_cent(Num_cells),this%flux_int(Num_cells+1))
                x0=spatial_discr_obj%init_point
                this%flux_int(1)=real_poly_1D(flux_coeffs,x0)
                select type (spatial_discr_obj)     
                type is (mesh_1D_Euler_homog_c)
                    do i=1,Num_cells
                        this%flux_cent(i)=real_poly_1D(flux_coeffs,x0+spatial_discr_obj%Delta_x*(2*i-1)/2d0)
                        this%flux_int(i+1)=real_poly_1D(flux_coeffs,x0+spatial_discr_obj%Delta_x*i)
                    end do
                    !this%flux_int(Num_cells+1)=real_poly_1D(flux_coeffs,x0+spatial_discr_obj%Delta_x*Num_cells)
                type is (mesh_1D_Euler_heterog_c)
                    !x=spatial_discr_obj%Delta_x(1)/2d0
                    !this%flux(1)=real_poly_1D(flux_coeffs,x)
                    x=x0
                    do i=1,Num_cells
                        this%flux_cent(i)=real_poly_1D(flux_coeffs,x+spatial_discr_obj%Delta_x(i)/2d0)
                        this%flux_int(i+1)=real_poly_1D(flux_coeffs,x+spatial_discr_obj%Delta_x(i))
                        x=x+spatial_discr_obj%Delta_x(i)
                    end do
                end select
                !Delta_x_n=spatial_discr_obj%get_Cell_size(Num_cells)
                !dqn_dx=der_real_poly_1D(flux_coeffs,spatial_discr_obj%measure-Delta_x_n/2d0)
                !d2qn_dx2=sec_der_real_poly_1D(flux_coeffs,spatial_discr_obj%measure-Delta_x_n/2d0)
                !if (this%source_term_order.eq.1) then
                !    q_out=this%flux(Num_cells) + Delta_x_n*dqn_dx/2d0 + (Delta_x_n**2)*d2qn_dx2/8d0 !> Taylor
                !end if
            !else if (spatial_discr_obj%scheme == 2) then
            !    allocate(this%flux(Num_cells+1))
            !    select type (spatial_discr_obj)     
            !    type is (mesh_1D_Euler_homog_c)
            !        do i=1,Num_cells+1
            !            this%flux(i)=real_poly_1D(flux_coeffs,spatial_discr_obj%Delta_x*(i-1))
            !        end do
            !    type is (mesh_1D_Euler_heterog_c)
            !        x=-spatial_discr_obj%Delta_x(1)
            !        do i=1,Num_cells+1
            !            x=x+spatial_discr_obj%Delta_x(i)
            !            this%flux(i)=real_poly_1D(flux_coeffs,x)
            !        end do
            !    end select
            !    q_out=this%flux(Num_cells+1)
            !else
            !    error stop "subroutine 'compute_flux_nonlin' not implemented yet for this scheme"
            !end if
        end subroutine
        
        subroutine compute_source_term(this,spatial_discr_obj,flux_coeffs) !> computes sink/source term from flux polynomial
            implicit none
            class(tpt_props_heterog_c) :: this
            class(spatial_discr_c), intent(in) :: spatial_discr_obj
            real(kind=8), intent(in) :: flux_coeffs(:) !> orden decreciente
            
            integer(kind=4) :: i
            real(kind=8) :: dq_dx
            allocate(this%source_term(spatial_discr_obj%Num_targets-spatial_discr_obj%targets_flag))
            
            select type (spatial_discr_obj)
            type is (mesh_1D_Euler_homog_c)
                if (spatial_discr_obj%scheme == 1 .and. spatial_discr_obj%targets_flag == 0) then
                    if (size(flux_coeffs) == 2) then
                        dq_dx=flux_coeffs(1) !> linear flux
                        do i=1,size(this%flux_cent)
                            this%source_term(i)=dq_dx
                        end do
                    else if (size(flux_coeffs) == 3) then
                        do i=1,size(this%flux_cent)
                            dq_dx=der_real_poly_1D(flux_coeffs,spatial_discr_obj%Delta_x*(2*i-1)/2d0) !> quadratic flux
                            this%source_term(i)=dq_dx
                        end do
                    else
                        error stop "Subroutine 'compute_source_term' not implemented yet for cubic fluxes"
                    end if
                else if (spatial_discr_obj%scheme == 2 .and. spatial_discr_obj%targets_flag == 0) then
                    do i=1,size(this%source_term)
                        this%source_term(i)=(this%flux_int(i+1)-this%flux_int(i))/spatial_discr_obj%Delta_x
                    end do
                else
                    error stop "Subroutine 'compute_source_term' not implemented yet for this scheme"
                end if
            end select
        end subroutine
        
        subroutine are_tpt_props_homog(this)
            implicit none
            class(tpt_props_heterog_c) :: this
            
            integer(kind=4) :: i
            real(kind=8), parameter :: eps=1d-12
            
            call are_diff_props_homog(this)
            if (this%homog_flag.eqv..true.) then
                do i=2,size(this%flux_int)
                    if (abs(this%flux_int(1)-this%flux_int(i))>eps) then
                        this%homog_flag=.false.
                        exit
                    end if
                end do
            end if
        end subroutine
        
        !subroutine allocate_flux_change(this,Num_chg)
        !implicit none
        !class(tpt_props_heterog_c) :: this !> heterogeneous transport properties object
        !integer(kind=4), intent(in) :: Num_chg !> number of time steps where flux changes
        !if (allocated(this%flux_change)) then
        !    error stop "Flux change already allocated"
        !end if
        !this%num_flux_chg=Num_chg
        !allocate(this%flux_change(this%num_flux_chg))
        !this%flux_change=0 !> initialize to zero
        !end subroutine
        
        !subroutine read_flux_change(this,filename)
        !implicit none
        !class(tpt_props_heterog_c) :: this !> heterogeneous transport properties object
        !character(len=*), intent(in) :: filename !> file name
        !
        !integer(kind=4) :: i,Num_chg
        !
        !open(unit=1,file=filename,status='old',action='read')
        !read(1,*) Num_chg
        !!if (allocated(this%flux_change)) then
        !!    error stop "Flux change already allocated"
        !!end if
        !call this%allocate_flux_change(Num_chg)
        !do i=1,this%num_flux_chg
        !    read(1,*) this%flux_change(i)
        !end do
        !close(1)
        !end subroutine
        
        subroutine update_flux_int_trans(this,k)
        implicit none
        class(tpt_props_heterog_c) :: this !> heterogeneous transport properties object
        integer(kind=4), intent(in) :: k !> new time step

        this%flux_int=this%flux_trans%time_series(k) !> update flux with the transient flux at time step k
        end subroutine
        
        subroutine set_cst_flux_flag(this,flag)
        implicit none
        class(tpt_props_heterog_c) :: this !> heterogeneous transport properties object
        logical, intent(in) :: flag !> TRUE if flux is constant, FALSE otherwise
        this%cst_flux_flag=flag !> set constant flux flag
        end subroutine
        
        

        subroutine compute_flux_rad(this,mesh,Q) !> computes radial flux
        implicit none
        class(tpt_props_heterog_c) :: this !> heterogeneous transport properties object
        class(spatial_discr_rad_c) :: mesh !> radial mesh object
        real(kind=8), intent(in) :: Q !> caudal
        !real(kind=8), intent(in) :: b !> thickness of the domain
        integer(kind=4) :: i !> loop index
        real(kind=8) :: r_prev,r_i !> radial coordinates
        real(kind=8), parameter :: pi=3.14159265358979323846 !> pi constant

        r_i=mesh%r_min !> initial radial coordinate
        !r_i=r_prev !> initial radial coordinate
        this%flux_int(1)=Q/(2d0*pi*mesh%r_min*mesh%targets(1)%thickness) !> radial flux at r_min
        do i=1,mesh%Num_targets !> loop over cells
            !r_i=r_prev+mesh%Delta_r(i)/2d0 !> radial coordinate
            this%flux_cent(i)=Q/(2d0*pi*(r_i+mesh%Delta_r(i)/2d0)*mesh%targets(i)%thickness) !> radial flux at interface
            r_i=r_i+mesh%Delta_r(i) !> update radial coordinate
            this%flux_int(i+1)=Q/(2d0*pi*r_i*mesh%targets(i)%thickness) !> radial flux at centre
        end do
        end subroutine

        subroutine compute_dispersion(this,scheme) !> computes dispersion from dispersivity and flux
            implicit none
            class(tpt_props_heterog_c) :: this !> heterogeneous transport properties object
            integer(kind=4), intent(in) :: scheme !> spatial discretization scheme
            !class(spatial_discr_c), intent(in) :: mesh !> spatial discretization object
            
            integer(kind=4) :: i
            
            !if (allocated(this%dispersivity) .and. allocated(this%flux_int)) then
            if (scheme==1) then !> traditional FD scheme
                if (allocated(this%dispersivity) .and. allocated(this%flux_cent)) then
                    allocate(this%dispersion(size(this%flux_cent)))
                    this%dispersion=this%dispersivity*this%flux_cent
                else if (allocated(this%dispersivity)) then
                    error stop "Flux not allocated"
                else
                    error stop "Dispersivity not allocated"
                end if
            else if (scheme==2) then !> proposed scheme
                if (allocated(this%dispersivity) .and. allocated(this%flux_int)) then
                    allocate(this%dispersion(size(this%flux_int)))
                    this%dispersion=this%dispersivity*this%flux_int
                else if (allocated(this%dispersivity)) then
                    error stop "Flux not allocated"
                else
                    error stop "Dispersivity not allocated"
                end if
            else
                error stop "Subroutine 'compute_dispersion' not implemented yet for this scheme"
            end if
        end subroutine
        
        subroutine allocate_flux(this,Num_nodes)
        class(tpt_props_heterog_c) :: this !> heterogeneous transport properties object
        integer(kind=4), intent(in) :: Num_nodes !> number of nodes in the mesh
        if (allocated(this%flux_int) .or. allocated(this%flux_cent)) then
            deallocate(this%flux_int,this%flux_cent)
        else
            allocate(this%flux_int(Num_nodes))
            allocate(this%flux_cent(Num_nodes-1))
        end if
        end subroutine
end module