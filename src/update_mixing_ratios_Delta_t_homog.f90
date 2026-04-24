!> Updates mixing ratios matrix with uniform time stepping
subroutine update_mixing_ratios_Delta_t_homog(this,Delta_t)
    use BCs_subroutines_m, only: Dirichlet_BCs_PDE, Dirichlet_Neumann_BCs_PDE, Neumann_homog_BCs, Robin_Neumann_homog_BCs, & 
        Robin_Dirichlet_BCs
    use transport_transient_m, only: transport_1D_transient_c
    use time_discr_m, only: time_discr_c, time_discr_homog_c
    use arrays_m, only: tridiag_matrix_c, prod_tridiag_mat_mat
    use metodos_sist_lin_m, only: compute_inverse_tridiag_matrix
    implicit none
    
    class(transport_1D_transient_c) :: this
    real(kind=8), intent(in) :: Delta_t
    
    integer(kind=4) :: i,j
    integer(kind=4) :: num_mix_rat !> number of mixing ratios in each target
    integer(kind=4) :: num_inf_bd_mix_rat !> number of inflow boundary mixing ratios outside the mesh
    integer(kind=4) :: num_out_bd_mix_rat !> number of outflow boundary mixing ratios outside the mesh
    integer(kind=4) :: num_rech_mix_rat !> number of recharge mixing ratios in each target
    real(kind=8) :: theta
    real(kind=8), parameter :: tol_inv=1d-12
    real(kind=8), allocatable :: A_mat_inv(:,:)
    
    type(tridiag_matrix_c) :: E_mat,X_mat_T,E_mat_prev
    
    !call this%mixing_ratios_conc%allocate_array(this%spatial_discr%Num_targets)
    !call this%mixing_ratios_R%allocate_array(this%spatial_discr%Num_targets)
    !call this%mixing_ratios_R_init%allocate_array(this%spatial_discr%Num_targets)
    
    num_mix_rat=0 !> we initialise total number of mixing ratios in each target
    !if (this%spatial_discr%targets_flag==0) then !> targets are cells
        if (this%BCs%labels(1).eq.1 .and. this%spatial_discr%targets_flag==0) then !> inflow Dirichlte BC, cell centred
            num_inf_bd_mix_rat=1
        !else if (this%BCs%labels(1).eq.1) then !> inflow Dirichlet BC, edge centred
        !    num_inf_bd_mix_rat=0
        else if (this%BCs%labels(1).eq.3) then !> inflow Robin BC
            num_inf_bd_mix_rat=1
        else
            num_inf_bd_mix_rat=0
        end if
        if (this%BCs%labels(2).eq.1 .and. this%spatial_discr%targets_flag==0) then !> outflow Dirichlet BC, cell centred
            num_out_bd_mix_rat=1
        else if (this%BCs%labels(2).eq.3) then !> outflow Robin BC
            num_out_bd_mix_rat=1
        else
            num_out_bd_mix_rat=0
        end if
    !else
    !    num_inf_bd_mix_rat=0
    !    num_out_bd_mix_rat=0
    !end if
    select type (this)
    class is (transport_1D_transient_c)
        if (this%tpt_props_heterog%cst_flux_flag .eqv. .false.) then
            num_rech_mix_rat=1 !> recharge term
        else
            num_rech_mix_rat=0
        end if
    end select
    num_mix_rat=num_mix_rat+num_inf_bd_mix_rat+num_out_bd_mix_rat+num_rech_mix_rat !> boundary and recharge mixing ratios
    
    ! select type (time_discr=>this%time_discr)
    ! class is (time_discr_homog_c)
    !     Delta_t=time_discr%Delta_t
    ! class default
    !     error stop "This subroutine is only applied if time discretisation is homogeneous"
    ! end select
    !Delta_t=this%time_discr%get_Delta_t()
!> We compute PDE arrays
    call this%compute_trans_mat_PDE()
    call this%compute_rech_mat_PDE()
    call this%compute_source_term_PDE()
    call this%compute_F_mat_PDE()
    !> We impose the boundary conditions
    if (this%BCs%labels(1).eq.1 .and. this%BCs%labels(2).eq.1) then
        call Dirichlet_BCs_PDE(this)
        !this%mixing_ratios_conc%cols(1)%dim=1+num_rech_mix_rat+num_inf_bd_mix_rat
        !this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%dim=1+num_rech_mix_rat+num_out_bd_mix_rat
        !num_Dir=2
    else if (this%BCs%labels(1).eq.1 .and. this%BCs%labels(2).eq.2) then
        call Dirichlet_Neumann_BCs_PDE(this)
        !this%mixing_ratios_conc%cols(1)%dim=1+num_rech_mix_rat+num_inf_bd_mix_rat
        !this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%dim=2+num_rech_mix_rat
    else if (this%BCs%labels(1).eq.2 .and. this%BCs%labels(2).eq.2) then
        call Neumann_homog_BCs(this)
        !this%mixing_ratios_conc%cols(1)%dim=2+num_rech_mix_rat
        !this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%dim=2+num_rech_mix_rat
    else if (this%BCs%labels(1).eq.3 .and. this%BCs%labels(2).eq.2) then
        call Robin_Neumann_homog_BCs(this)
        !this%mixing_ratios_conc%cols(1)%dim=2+num_rech_mix_rat
        !this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%dim=2+num_rech_mix_rat
    else if (this%BCs%labels(1).eq.3 .and. this%BCs%labels(2).eq.1) then
        call Robin_Dirichlet_BCs(this)
        !this%mixing_ratios_conc%cols(1)%dim=2+num_rech_mix_rat
        !this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%dim=1+num_rech_mix_rat+num_out_bd_mix_rat
    else
        error stop "Boundary conditions not implemented yet"
    end if
    if (this%time_discr%int_method.eq.1) then !> Euler explicit
        theta=0d0
        !do i=2,this%mixing_ratios_conc%num_cols-1
        !    this%mixing_ratios_conc%cols(i)%dim=3+num_rech_mix_rat
        !    this%mixing_ratios_R%cols(i)%dim=1
        !    this%mixing_ratios_R_init%cols(i)%dim=1
        !end do
        !call this%mixing_ratios_conc%allocate_columns()
        !call this%mixing_ratios_R%allocate_columns()
        !call this%mixing_ratios_R_init%allocate_columns()
        !call this%allocate_mix_conc_indices()
    else !> Euler implicit
        num_mix_rat=num_mix_rat+this%spatial_discr%Num_targets
        !call this%mix_conc_indices%allocate_array(this%mixing_ratios_conc%num_cols)
        !do i=1,this%mixing_ratios_conc%num_cols
        !    call this%mixing_ratios_conc%cols(i)%set_dim(num_mix_rat)
        !    call this%mixing_ratios_conc%cols(i)%allocate_vector()
        !    call this%mixing_ratios_R%cols(i)%set_dim(this%mixing_ratios_R%num_cols)
        !    call this%mixing_ratios_R_init%cols(i)%set_dim(this%mixing_ratios_R%num_cols)
        !    call this%mixing_ratios_R%cols(i)%allocate_vector()
        !    call this%mixing_ratios_R_init%cols(i)%allocate_vector()
        !    !call this%mix_conc_indices%cols(i)%allocate_vector(this%mixing_ratios_conc%cols(i)%dim+2)
        !end do
        if (this%time_discr%int_method.eq.2) then !> Euler fully implicit
            theta=1d0
        else if (this%time_discr%int_method.eq.3) then !> Crank-Nicolson
            theta=5d-1
        else
            error stop "Time discretisation not implemented yet"
        end if
    end if
    !call this%allocate_mix_conc_indices()
!> We compute arrays for linear system
    call this%compute_E_mat_1D(Delta_t,E_mat,E_mat_prev)
    call this%compute_X_mat_1D(theta,E_mat_prev)
    call this%compute_A_mat_1D(theta,E_mat)
    call this%compute_Y_mat()
    call this%compute_Z_mat()
    call this%compute_f_vec(Delta_t)
    
    if (theta>0d0) then !> Euler implicit
        !call this%allocate_A_mat_inv() !> allocate A_mat_inv
        !call this%allocate_mixing_ratios_mat_conc_mesh()
        !allocate(this%mixing_ratios_mat_conc(this%mixing_ratios_conc%num_cols,this%mixing_ratios_conc%num_cols))
        call compute_inverse_tridiag_matrix(this%A_mat,tol_inv,this%A_mat_inv)
        if (theta.eq.1d0) then !> Euler fully implicit
            this%mixing_ratios_mat_conc_mesh=transpose(this%A_mat_inv)
        else
            call this%X_mat%compute_transpose_tridiag_matrix(X_mat_T)
            this%mixing_ratios_mat_conc_mesh=prod_tridiag_mat_mat(X_mat_T,transpose(this%A_mat_inv))
        end if
        !this%mixing_ratios_mat_Rk=transpose(A_mat_inv) !> mixing ratios of kinetic reaction amounts
        !if (num_mix_rat>2*this%spatial_discr%Num_targets .or. num_mix_rat<2*this%spatial_discr%num_targets) then !> boundary and no recharge
        !    this%mixing_ratios_conc%cols(1)%col_1=this%bd_mat(1)*this%A_mat_inv(:,1) !> inflow boundary mixing ratios
        !    this%mixing_ratios_conc%cols(this%spatial_discr%Num_targets+2)%col_1=this%bd_mat(2)*this%A_mat_inv(:,this%spatial_discr%num_targets) !> outflow boundary mixing ratios
        !else if (num_mix_rat>this%spatial_discr%Num_targets) then !> recharge with no boundary
        !    this%mixing_ratios_mat_Rk=this%A_mat_inv
        !end if
        !if (num_inf_bd_mix_rat==1) then !> inflow boundary mixing ratios
        !    this%mixing_ratios_conc%cols(1)%col_1=this%bd_mat(1)*this%A_mat_inv(:,1)
        !    this%mix_conc_indices%cols(1)%col_1(1)=1 !> current water
        !    this%mix_conc_indices%cols(1)%col_1(this%mix_conc_indices%cols(1)%dim-1)=0 !> upstream waters
        !    this%mix_conc_indices%cols(1)%col_1(this%mix_conc_indices%cols(1)%dim)=this%spatial_discr%num_targets+num_out_bd_mix_rat !> downstream waters
        !end if
        !if (num_out_bd_mix_rat==1) then !> outflow boundary mixing ratios
        !    this%mixing_ratios_conc%cols(this%spatial_discr%num_targets+num_inf_bd_mix_rat+1)%col_1=this%bd_mat(2)*this%A_mat_inv(:,this%spatial_discr%Num_targets)
        !    this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_bd_mix_rat+1)%col_1(1)=this%spatial_discr%num_targets+num_inf_bd_mix_rat+1 !> current water
        !    this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_bd_mix_rat+1)%col_1(this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_bd_mix_rat+1)%dim-1)=this%spatial_discr%num_targets+num_inf_bd_mix_rat !> upstream waters
        !    this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_bd_mix_rat+1)%col_1(this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_bd_mix_rat+1)%dim)=0 !> downstream waters
        !end if
        do i=1,this%spatial_discr%num_targets
            this%mixing_ratios_conc%cols(i)%col_1(1+num_inf_bd_mix_rat)=this%Z_mat(1)*this%A_mat_inv(i,1) !> concentration mixing ratio of inflow boundary
            this%mixing_ratios_conc%cols(i)%col_1(this%spatial_discr%num_targets+num_inf_bd_mix_rat+num_out_bd_mix_rat)=&
                this%Z_mat(2)*this%A_mat_inv(i,this%spatial_discr%num_targets) !> concentration mixing ratio of outflow boundary
            this%mixing_ratios_conc%cols(i)%col_1(1)=this%mixing_ratios_mat_conc_mesh(i,i) !> concentration mixing ratio of current water
            this%mixing_ratios_R%cols(i)%col_1(1)=this%A_mat_inv(i,i) !> reaction amount mixing ratio of current water
            !this%mixing_ratios_R_init%cols(i)%col_1(1)=this%mixing_ratios_R%cols(i)%col_1(1) !> initial reaction amount mixing ratio of current water
            this%mix_conc_indices%cols(i)%col_1(1+num_inf_bd_mix_rat)=num_inf_bd_mix_rat !> inflow water index
            this%mix_conc_indices%cols(i)%col_1(this%spatial_discr%num_targets+num_inf_bd_mix_rat+num_out_bd_mix_rat)=&
                this%spatial_discr%num_targets+num_inf_bd_mix_rat+num_out_bd_mix_rat !> outflow water index
            this%mix_conc_indices%cols(i)%col_1(1)=num_inf_bd_mix_rat+i !> current water index
            this%mix_react_indices%cols(i)%col_1(1)=this%mix_conc_indices%cols(i)%col_1(1) !> current water index
            this%mix_conc_indices%cols(i)%col_1(this%mix_conc_indices%cols(i)%dim-1)=num_inf_bd_mix_rat+i-1 !> upstream waters
            this%mix_react_indices%cols(i)%col_1(this%mix_react_indices%cols(i)%dim-1)=i-1 !> upstream waters
            this%mix_conc_indices%cols(i)%col_1(this%mix_conc_indices%cols(i)%dim)=this%mixing_ratios_conc%cols(i)%dim-&
                i-num_inf_bd_mix_rat !> downstream waters
            this%mix_react_indices%cols(i)%col_1(this%mix_react_indices%cols(i)%dim)=&
                this%spatial_discr%num_targets-i !> downstream waters
            do j=1,i-1
                this%mixing_ratios_conc%cols(i)%col_1(1+num_inf_bd_mix_rat+j)=this%mixing_ratios_mat_conc_mesh(j,i)
                this%mixing_ratios_R%cols(i)%col_1(1+j)=this%A_mat_inv(i,j)
                !this%mixing_ratios_R_init%cols(i)%col_1(1+j)=this%mixing_ratios_R%cols(i)%col_1(1+j) !> initial reaction amount mixing ratio of downstream waters
                this%mix_conc_indices%cols(i)%col_1(1+num_inf_bd_mix_rat+j)=num_inf_bd_mix_rat+j
                this%mix_react_indices%cols(i)%col_1(1+j)=this%mix_conc_indices%cols(i)%col_1(1+num_inf_bd_mix_rat+j)
            end do
            do j=i+1,this%spatial_discr%num_targets
                this%mixing_ratios_conc%cols(i)%col_1(num_inf_bd_mix_rat+j)=this%mixing_ratios_mat_conc_mesh(j,i)
                this%mixing_ratios_R%cols(i)%col_1(j)=this%A_mat_inv(i,j)
                !this%mixing_ratios_R_init%cols(i)%col_1(j)=this%mixing_ratios_R%cols(i)%col_1(j) !> initial reaction amount mixing ratio of upstream waters
                this%mix_conc_indices%cols(i)%col_1(num_inf_bd_mix_rat+j)=num_inf_bd_mix_rat+j
                this%mix_react_indices%cols(i)%col_1(j)=this%mix_conc_indices%cols(i)%col_1(num_inf_bd_mix_rat+j)
            end do
            do j=1,num_rech_mix_rat
                this%mixing_ratios_conc%cols(i)%col_1(this%spatial_discr%num_targets+num_inf_bd_mix_rat+num_out_bd_mix_rat+j)=&
                    this%Y_mat%diag(j)*this%A_mat_inv(i,j) !> mixing ratio of recharge waters
                this%mix_conc_indices%cols(i)%col_1(this%spatial_discr%num_targets+num_inf_bd_mix_rat+num_out_bd_mix_rat+j)=&
                    this%spatial_discr%num_targets+num_inf_bd_mix_rat+num_out_bd_mix_rat+j !> recharge water index
            end do
        end do
    else !> Euler explicit
        !> Mixing ratios of reaction amounts is the identity matrix
        this%mix_conc_indices%cols(1)%col_1(1)=1+num_inf_bd_mix_rat
        if (this%BCs%labels(1).eq.1) then !> Dirichlet BC inflow
            this%mixing_ratios_conc%cols(1)%col_1(1+num_inf_bd_mix_rat+num_rech_mix_rat)=this%Y_mat%diag(1)
            this%mixing_ratios_conc%cols(1)%col_1(1+num_inf_bd_mix_rat)=this%Z_mat(1)
            this%mixing_ratios_conc%cols(1)%col_1(1)=this%X_mat%diag(1) 
            this%mix_conc_indices%cols(1)%col_1(2)=num_inf_bd_mix_rat
            this%mix_conc_indices%cols(1)%col_1(3)=this%spatial_discr%num_targets+num_inf_bd_mix_rat+num_out_bd_mix_rat+&
                num_rech_mix_rat
        else 
            this%mixing_ratios_conc%cols(1)%col_1=[this%X_mat%diag(1),this%X_mat%super(1)] !> mixing ratios in current water
            this%mix_conc_indices%cols(1)%col_1(2)=2 !> downstream water
            this%mix_conc_indices%cols(1)%col_1(3)=0 !> number of upstream waters
            this%mix_conc_indices%cols(1)%col_1(4)=1 !> number of downstream waters
        end if
        this%mixing_ratios_R%cols(1)%col_1=1d0 !> mixing ratio of kinetic reaction rate in first target
        !this%mixing_ratios_R_init%cols(1)%col_1=this%mixing_ratios_R%cols(1)%col_1 !> initial mixing ratio of kinetic reaction rate in first target
        do i=2,this%mixing_ratios_conc%num_cols-1
            this%mixing_ratios_conc%cols(i)%col_1=[this%X_mat%diag(i),this%X_mat%sub(i-1),this%X_mat%super(i)]
            this%mixing_ratios_R%cols(i)%col_1=1d0 !> mixing ratio of kinetic reaction rate in i-th target
            !this%mixing_ratios_R_init%cols(i)%col_1=this%mixing_ratios_R%cols(i)%col_1 !> initial mixing ratio of kinetic reaction rate in i-th target
            this%mix_conc_indices%cols(i)%col_1(1)=i !> current water
            this%mix_conc_indices%cols(i)%col_1(2)=i-1 !> upstream water
            this%mix_conc_indices%cols(i)%col_1(3)=i+1 !> downstream water
            this%mix_conc_indices%cols(i)%col_1(4)=1 !> number of upstream waters
            this%mix_conc_indices%cols(i)%col_1(5)=1 !> number of downstream waters
        end do
        this%mixing_ratios_R%cols(this%mixing_ratios_conc%num_cols)%col_1=1d0 !> mixing ratio of kinetic reaction rate in last target
        !this%mixing_ratios_R_init%cols(this%mixing_ratios_conc%num_cols)%col_1=&
        !    this%mixing_ratios_R%cols(this%mixing_ratios_conc%num_cols)%col_1 !> initial mixing ratio of kinetic reaction rate in last target
        this%mix_conc_indices%cols(this%mixing_ratios_conc%num_cols)%col_1(1)=1
        if (this%BCs%labels(2).eq.1) then !> Dirichlet BC outflow
            this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%col_1(this%mixing_ratios_conc%cols(&
                this%mixing_ratios_conc%num_cols)%dim)=&
                this%X_mat%diag(this%mixing_ratios_conc%num_cols) !> mixing ratio of current water
            this%mix_conc_indices%cols(this%mixing_ratios_conc%num_cols)%col_1(2)=0 !> number of upstream waters
            this%mix_conc_indices%cols(1)%col_1(3)=0 !> number of downstream waters
        else
            this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%col_1=[&
                this%X_mat%diag(this%mixing_ratios_conc%num_cols),this%X_mat%sub(this%mixing_ratios_conc%num_cols-1)]
            this%mix_conc_indices%cols(this%mixing_ratios_conc%num_cols)%col_1(2)=this%mixing_ratios_conc%num_cols-1 !> upstream water
            this%mix_conc_indices%cols(this%mixing_ratios_conc%num_cols)%col_1(3)=1 !> number of upstream waters
            this%mix_conc_indices%cols(this%mixing_ratios_conc%num_cols)%col_1(4)=0 !> number of downstream waters
        end if
    end if
end subroutine 