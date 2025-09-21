!> Writes transport data from 1D reactive transport model
!> File is already opened
subroutine write_transport_data(this,unit)
    !use analytical_solutions_transport_m
    use RT_1D_m, only: RT_1D_c, RT_1D_transient_c
    use time_discr_m, only: time_discr_c
    !use chem_system_m, only: chem_system_c
    use spatial_discr_rad_m
    implicit none
    class(RT_1D_c), intent(in) :: this
    integer(kind=4), intent(in) :: unit
    !character(len=*), intent(in) :: file_out
    
    class(time_discr_c), pointer :: time_discr
    integer(kind=4) :: i,j,n,num_cells
    integer(kind=4), allocatable :: inz(:,:)
    real(kind=8), allocatable :: props_mat(:,:),Delta_t(:),stab_params(:),anal_sol(:,:),gamma_1(:),gamma_2(:),c2(:,:)
    !class(chem_system_c), pointer :: my_chem_syst
    
    
    
    !open(unit=unit,file=file_out,status='unknown',form='formatted')
    
    select type (this)
    type is (RT_1D_transient_c)
        !write(unit,"(2x,'Equation:',5x,'U*dc/dt*F = U*c_nc + U*S_nc,k^T*r_k',/)")
        n=this%transport%spatial_discr%Num_targets
        num_cells=this%transport%spatial_discr%Num_targets-this%transport%spatial_discr%targets_flag !> number of cells
        if (this%int_method_chem_reacts == 1) then
            write(unit,"(2x,'Integration method chemical reactions:',10x,'Euler explicit',/)")
        !else if (this%transport%time_discr%int_method == 2) then
        !>    write(unit,"(2x,'Integration method:',10x,'Euler semi-implicit',/)")
        else if (this%int_method_chem_reacts == 2) then
            write(unit,"(2x,'Integration method chemical reactions:',10x,'Euler fully implicit',/)")
        !else if (this%transport%time_discr%int_method.eqv.4) then
        !>    write(unit,"(2x,'Integration method chemical reactions:',10x,'Crank-Nicolson',/)")
        else
            error stop "Integration method not implemented yet for RT"
        end if
        !select type (mesh=>this%transport%spatial_discr)
        !type is (spatial_discr_rad_c)
        !    write(unit,"(2x,'Dimension:',I5/)") mesh%dim
        !    if (mesh%dim.eqv.1) then
        !        write(unit,"(2x,'Length of domain:',F15.5/)") mesh%measure
        !    else
        !        write(unit,"(2x,'Radius:',F15.5/)") mesh%radius
        !    end if
        !type is (mesh_1D_Euler_homog_c)
        !    write(unit,"(2x,'Length of domain:',F15.5/)") mesh%measure
        !end select
        !write(unit,"(2x,'Number of cells:',I5/)") num_cells
        write(unit,"(2x,'Number of targets:',I5/)") this%transport%spatial_discr%Num_targets!-this%spatial_discr%targets_flag
        !if (this%transport%spatial_discr%targets_flag.eqv.0) then
        !    write(unit,"(2x,'Targets type:',10x,'Cells',/)")
        !else
        !    write(unit,"(2x,'Targets type:',10x,'Interfaces',/)")
        !end if
        !if (this%transport%BCs%BCs_label(1).eqv.1 .and. this%transport%BCs%BCs_label(2).eqv.1) then
        !    write(unit,"(2x,'Boundary conditions:',10x,'Dirichlet',/)")
        !    !write(unit,"(2x,'Boundary conditions:',10x,'Dirichlet',5x,ES15.5,5x,ES15.5/)") this%transport%BCs%conc_inf, this%transport%BCs%conc_out
        !else if (this%transport%BCs%BCs_label(1).eqv.2 .and. this%transport%BCs%BCs_label(2).eqv.2) then
        !    write(unit,"(2x,'Boundary conditions:',10x,'Neumann homogeneous',/)")
        !else if (this%transport%BCs%BCs_label(1).eqv.3 .and. this%transport%BCs%BCs_label(2).eqv.2) then
        !    write(unit,"(2x,'Boundary conditions:',10x,'Prescribed mass flux at inflow',/)")
        !    !write(unit,"(2x,'Boundary conditions:',10x,'Prescribed mass flux',5x,'Flux inflow:',ES15.5,5x,'Concentration inflow:',ES15.5/)") this%transport%BCs%flux_inf, this%transport%BCs%conc_inf
        !end if
        !if (this%transport%spatial_discr%scheme.eqv.1) then
        !>    write(unit,"(2x,'Scheme:',10x,'CFD',/)")
        !else if (this%transport%spatial_discr%scheme.eqv.2) then
        !>    write(unit,"(2x,'Scheme:',10x,'IFD',/)")
        !else if (this%transport%spatial_discr%scheme.eqv.3) then
        !>    write(unit,"(2x,'Scheme:',10x,'Upwind',/)")
        !else
        !>    error stop "Scheme not implemented yet"
        !end if
        !select type (transport=>this%transport)
        !type is (transport_1D_transient_c)
            !write(unit,"(2x,'Transport properties:'/)")
            !write(unit,"(10x,'Porosity:',2x,*(ES15.5)/)") this%transport%tpt_props_heterog%porosity
            !write(unit,"(10x,'Dispersion:',2x,*(ES15.5)/)") this%transport%tpt_props_heterog%dispersion
            !write(unit,"(10x,'Flux:',2x,*(ES15.5)/)") this%transport%tpt_props_heterog%flux
            write(unit,"(/,2x,'Time step:'/)")
            write(unit,"(2x,ES15.5/)") this%transport%time_discr%get_Delta_t()
            write(unit,"(2x,'Final time:'/)")
            write(unit,"(2x,ES15.5/)") this%transport%time_discr%Final_time
            !write(unit,"(/,2x,'Stability parameters:'/)")
            !write(unit,"(10x,'Critical time step for transport:',ES15.5,/)") this%transport%stab_params_tpt%Delta_t_crit
            !write(unit,"(10x,'Courant:',ES15.5/)") this%transport%stab_params_tpt%Courant
            !write(unit,"(10x,'Peclet:',ES15.5,/)") this%transport%stab_params_tpt%Peclet
            write(unit,"(/,2x,'Dimension + Mixing ratios (by rows) (including boundary and sink/source terms):'/)")
            do i=1,this%transport%mixing_ratios_conc%num_cols
                write(unit,"(2x,I5,*(ES15.5))") this%transport%mixing_ratios_conc%cols(i)%dim, (&
                this%transport%mixing_ratios_conc%cols(i)%col_1(j), j=1,this%transport%mixing_ratios_conc%cols(i)%dim)
            end do
            write(unit,"(/,2x,'Mixing waters indices:'/)")
            do i=1,this%transport%mixing_waters_indices%num_cols
                write(unit,"(2x,I5,*(I5))") this%chemistry%num_ext_waters+i, this%transport%mixing_waters_indices%cols(i)%col_1
            end do
            !if (this%transport%time_discr%int_method.eqv.1) then
            !    !write(unit,"(10x,'The first column contains the diagonal elements in the mixing ratios matrix'/)")
            !    do i=1,this%transport%mixing_ratios_conc%num_cols
            !        write(unit,"(2x,*(F15.5))") (this%transport%mixing_ratios_conc%cols(i)%col_1(j), j=1,this%transport%mixing_ratios_conc%cols(i)%dim)
            !    end do
            !else if (this%transport%time_discr%int_method.eqv.2) then
            !    !write(unit,"(2x,*(F15.5))") (this%transport%mixing_ratios_conc%cols(1)%col_1(j), j=1,this%transport%mixing_ratios_conc%num_cols)
            !    !do i=2,this%transport%mixing_ratios_conc%num_cols-1
            !    !    write(unit,"(2x,*(F15.5))") (this%transport%mixing_ratios_conc%cols(i)%col_1(j), j=2,i),  this%transport%mixing_ratios_conc%cols(i)%col_1(1), (this%transport%mixing_ratios%cols(i)%col_1(j), j=i+1,this%transport%mixing_ratios%num_cols)
            !    !end do
            !    !write(unit,"(2x,*(F15.5))") (this%transport%mixing_ratios%cols(this%transport%mixing_ratios%num_cols)%col_1(j), j=2,this%transport%mixing_ratios%num_cols), this%transport%mixing_ratios%cols(this%transport%mixing_ratios%num_cols)%col_1(1)
            !    do i=1,this%transport%mixing_ratios%num_cols
            !        write(unit,"(2x,*(F15.5))") (this%transport%mixing_ratios_mat(j,i), j=1,this%transport%mixing_ratios%num_cols)
            !    end do
            !end if
            !if (n.eqv.1) then
            !>    write(unit,"(2x,F15.5,/)") this%transport%mixing_ratios%(n,2), this%transport%mixing_ratios(n,4)
            !else if (n.eqv.2) then
            !>    write(unit,"(2x,F15.5,/)") (this%transport%mixing_ratios(1,j), j=2,4)
            !>    write(unit,"(2x,F15.5,/)") (this%transport%mixing_ratios(n,j),j=1,2), this%transport%mixing_ratios(n,4)
            !else
            !>    write(unit,"(2x,15x,3F15.5)") (this%transport%mixing_ratios(1,j), j=2,4)
            !>    do i=2,n-1
            !>        write(unit,"(2x,4F15.5)") (this%transport%mixing_ratios(i,j), j=1,4)!this%transport%mixing_ratios%sub(i-1), this%transport%mixing_ratios%diag(i), this%transport%mixing_ratios%super(i), this%transport%f_vec(i)
            !>    end do
            !>    write(unit,"(2x,2F15.5,15x,F15.5/)") (this%transport%mixing_ratios(n,j), j=1,2), this%transport%mixing_ratios(n,4)!this%transport%mixing_ratios%sub(n-1), this%transport%mixing_ratios%diag(n), this%transport%f_vec(n)
            !end if
            !if (this%transport%time_discr%int_method.eqv.1 .or. this%transport%time_discr%int_method.eqv.2) then
            !>    if (n.eqv.1) then
            !>        write(unit,"(2x,F15.5,/)") this%transport%X_mat%diag(1)
            !>    else if (n>2) then
            !>        write(unit,"(2x,15x,3F15.5)") this%transport%X_mat%diag(1), this%transport%X_mat%super(1), this%transport%f_vec(1)
            !>        do i=2,n-1
            !>            write(unit,"(2x,4F15.5)") this%transport%X_mat%sub(i-1), this%transport%X_mat%diag(i), this%transport%X_mat%super(i), this%transport%f_vec(i)
            !>        end do
            !>        write(unit,"(2x,2F15.5,15x,F15.5/)") this%transport%X_mat%sub(n-1), this%transport%X_mat%diag(n), this%transport%f_vec(n)
            !>    end if
            !else if (this%transport%time_discr%int_method.eqv.3) then
                !write(unit,"(/,2x,'Matrix that is inverted:'/)")
                !if (n.eqv.1) then
                !    write(unit,"(17x,ES15.5,/)") this%transport%A_mat%diag(1)
                !else if (n>2) then
                !    write(unit,"(17x,*(ES15.5))") this%transport%A_mat%diag(1), this%transport%A_mat%super(1)
                !    do i=2,n-1
                !        write(unit,"(2x,*(ES15.5))") this%transport%A_mat%sub(i-1), this%transport%A_mat%diag(i), this%transport%A_mat%super(i)
                !    end do
                !    write(unit,"(2x,*(ES15.5)/)") this%transport%A_mat%sub(n-1), this%transport%A_mat%diag(n)
                !end if
            !end if
    end select
    !close(unit)
end subroutine