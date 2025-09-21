!> Writes data and results of 1D transient transport equation
subroutine write_transport_1D_transient(this,root,Time_out,output)
    use transport_transient_m, only: transport_1D_transient_c
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c
    use time_discr_m, only: time_discr_homog_c
    implicit none
    !> Variables
    class(transport_1D_transient_c), intent(in) :: this !> 1D transient transport object
    character(len=*), intent(in) :: root !> root name for output files
    real(kind=8), intent(in) :: Time_out(:)
    real(kind=8), intent(in) :: output(:,:)

    integer(kind=4) :: Num_output
    real(kind=8), allocatable :: stab_params(:)
    integer(kind=4) :: i,j,k,n,n_flux
    character(len=256) :: file_out
    real(kind=8) :: Delta_t=0d0
    
    n=this%spatial_discr%Num_targets
    Num_output=size(Time_out)
    !if (this%spatial_discr%scheme.eq.1) then
    !    if (this%dimless.eqv..true.) then
    !        write(file_out,root//'transport_1D_trans_CFDS_adim.out')
    !    else
    !        write(file_out,"('transport_1D_trans_CFDS.out')")
    !    end if
    !else if (this%spatial_discr%scheme.eq.2) then
    !    if (this%dimless.eqv..true.) then
    !        write(file_out,"('transport_1D_trans_IFDS_adim.out')")
    !    else
    !        write(file_out,root//'transport_1D_trans_IFDS.out')
    !    end if
    !end if
    open(unit=1,file=root//'.out',status='unknown')
    write(1,"(2x,'Equation:',5x,'F*dc/dt = T*c + g',/)")
    write(1,"(2x,'Mesh volume:',F15.5/)") this%spatial_discr%measure
    write(1,"(2x,'Number of targets:',I5/)") this%spatial_discr%Num_targets!-this%spatial_discr%targets_flag
    if (this%spatial_discr%targets_flag.eq.0) then
        write(1,"(2x,'Targets type: Cells',/)")
    else
        write(1,"(2x,'Targets type: Interfaces',/)")
    end if
    select type (mesh=>this%spatial_discr)
    type is (mesh_1D_Euler_homog_c)
        write(1,"(2x,'Mesh size:',F15.5/)") mesh%Delta_x
    end select
    if (this%BCs%BCs_label(1).eq.1 .and. this%BCs%BCs_label(2).eq.1) then
        write(1,"(2x,'Boundary conditions:',10x,'Dirichlet',10x,2F15.5,/)") this%BCs%conc_inf, this%BCs%conc_out
    else if (this%BCs%BCs_label(1).eq.2 .and. this%BCs%BCs_label(2).eq.2) then
        write(1,"(2x,'Boundary conditions:',10x,'Neumann homogeneous',/)")
    else if (this%BCs%BCs_label(1).eq.3 .and. this%BCs%BCs_label(2).eq.2) then
        write(1,"(2x,'Boundary conditions:',10x,'Robin inflow, Neumann homogeneous outflow',/)")
    end if
    if (this%spatial_discr%scheme.eq.1) then
        write(1,"(2x,'Scheme:',10x,'Traditional FD',/)")
    else if (this%spatial_discr%scheme.eq.2) then
        write(1,"(2x,'Scheme:',10x,'Petchamé-Carrera (2024)',/)")
    else if (this%spatial_discr%scheme.eq.3) then
        write(1,"(2x,'Scheme:',10x,'Upwind',/)")
    end if
    select type (time=>this%time_discr)
    type is (time_discr_homog_c)
        Delta_t=this%time_discr%get_Delta_t()
        write(1,"(2x,'Time step:',ES15.5/)") Delta_t
    end select
    write(1,"(2x,'Final time:',ES15.5/)") this%time_discr%Final_time
    write(1,"(2x,'Number of time steps:',I10/)") this%time_discr%Num_time
    if (this%time_discr%int_method.eq.1) then
        write(1,"(2x,'Integration method:',10x,'Euler explicit',/)")
    else if (this%time_discr%int_method.eq.2) then
        write(1,"(2x,'Integration method:',10x,'Euler fully implicit',/)")
    else if (this%time_discr%int_method.eq.3) then
        write(1,"(2x,'Integration method:',10x,'Crank-Nicolson',/)")
    else if (this%time_discr%int_method.eq.4) then
        write(1,"(2x,'Integration method:',10x,'RKF45',/)")
    end if
    write(1,"(2x,'Linear system:',5x,'A*c^(k+1) = B*c^k + f',/)")
    if (this%dimless.eqv..false.) then
        write(1,"(2x,'Properties:'/)")
        write(1,"(10x,'Porosity:',/)")
        do i=1,n-this%spatial_discr%targets_flag
            write(1,"(20x,ES15.5)") this%tpt_props_heterog%porosity(i)
        end do
        write(1,"(/,10x,'Flux at interfaces:'/)")
        do i=1,size(this%tpt_props_heterog%flux_int)
            write(1,"(20x,ES15.5)") this%tpt_props_heterog%flux_int(i)
        end do
        write(1,"(10x,'Dispersion:',/)")
        do i=1,size(this%tpt_props_heterog%dispersion)
            write(1,"(20x,ES15.5)") this%tpt_props_heterog%dispersion(i)
        end do
    end if
    if (this%tpt_props_heterog%homog_flag.eqv..true.) then
        write(1,"(/,2x,'Stability parameters:'/)")
        write(1,"(10x,'Critical time step:',/)")
        write(1,"(20x,ES15.5)") this%stab_params_tpt%Delta_t_crit
        if (Delta_t>0d0 .and. Delta_t>this%stab_params_tpt%Delta_t_crit .and. this%time_discr%int_method<3) then
            write(1,"(/,20x,'You must reduce time step to have stability')")
        end if
        write(1,"(/,10x,'Courant:'/)")
        write(1,"(20x,ES15.5)") this%stab_params_tpt%Courant
        if (this%stab_params_tpt%Courant>1d0) then
            write(1,"(/,20x,'Courant condition violated'/)") 
        end if
        write(1,"(/,10x,'Peclet:'/)")
        write(1,"(20x,ES15.5)") this%stab_params_tpt%Peclet
        if (this%stab_params_tpt%Peclet>2d0) then
            write(1,"(/,20x,'Peclet condition violated'/)") 
        end if
    end if
    !write(1,"(2x,'F:'/)") 
    !do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
    !>    write(1,"(2x,F15.5)") this%F_mat%diag(i)
    !end do
    !write(1,"(/,2x,'Transition matrix T (with BCs):'/)") 
    !write(1,"(17x,2F15.5)") this%trans_mat%diag(1), this%trans_mat%super(1)    
    !do i=2,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag-1
    !>    write(1,"(2x,3F15.5)") this%trans_mat%sub(i-1), this%trans_mat%diag(i), this%trans_mat%super(i)
    !end do
    !write(1,"(2x,2F15.5/)") this%trans_mat%sub(this%spatial_discr%Num_targets-1), this%trans_mat%diag(this%spatial_discr%Num_targets)
    write(1,"(/,2x,'Matrix A:'/)") 
    write(1,"(17x,2F15.5)") this%A_mat%diag(1), this%A_mat%super(1)    
    do i=2,this%spatial_discr%Num_targets-1
        write(1,"(2x,3F15.5)") this%A_mat%sub(i-1), this%A_mat%diag(i), this%A_mat%super(i)
    end do
    write(1,"(2x,2F15.5/)") this%A_mat%sub(this%spatial_discr%Num_targets-1), this%A_mat%diag(this%spatial_discr%Num_targets)
    write(1,"(/,2x,'Matrix B:'/)") 
    write(1,"(17x,2F15.5)") this%X_mat%diag(1), this%X_mat%super(1)    
    do i=2,this%spatial_discr%Num_targets-1
        write(1,"(2x,3F15.5)") this%X_mat%sub(i-1), this%X_mat%diag(i), this%X_mat%super(i)
    end do
    write(1,"(2x,2F15.5/)") this%X_mat%sub(this%spatial_discr%Num_targets-1), this%X_mat%diag(this%spatial_discr%Num_targets)
    write(1,"(/,2x,'Vector f:'/)")
    do i=1,this%spatial_discr%Num_targets
        write(1,"(2x,F15.5)") this%f_vec(i)
    end do
   
    if (size(output,1).eq.this%spatial_discr%Num_targets-this%spatial_discr%targets_flag) then
        if (this%time_discr%int_method<4) then
            write(1,"(/,2x,'Cell',*(ES20.5)/)") (Time_out(k), k=1,Num_output)
            do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
                write(1,"(2x,I4,*(F20.5))") i,(output(i,k), k=1,Num_output)
            end do
        else
            write(1,"(/,2x,'Cell',2ES20.5/)") Time_out(1), this%time_discr%Final_time
            do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
                write(1,"(2x,I5,2ES20.10)") i, this%conc_init(i), this%conc_obj%conc(i)
            end do
        end if
    else
        write(1,"(/,2x,'Mobile zone:'/)")
        write(1,"(10x,'Cell',3ES20.5/)") (Time_out(k), k=1,Num_output)
        do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
            write(1,"(10x,I4,3F20.5)") i,(output(i,k), k=1,Num_output)
        end do
        write(1,"(/,2x,'Immobile zones:'/)")
        do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
           !> write(1,"(10x,I4,3F20.5)") i,(output(this%spatial_discr%Num_targets-this%spatial_discr%targets_flag+i,k), k=1,Num_output)
        end do
    end if
    rewind(1)
    close(1)
end subroutine write_transport_1D_transient