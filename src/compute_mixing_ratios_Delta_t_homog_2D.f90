!> \file compute_mixing_ratios_Delta_t_homog_2D.f90
!> \brief Computes mixing ratios matrix for reactive transport with uniform time stepping
!> \details This subroutine computes the mixing ratios that describe how target waters mix
!> with boundary waters, neighboring waters, and recharge waters during one time step.
!> For implicit methods, it computes and inverts the system matrix A to obtain mixing ratios.
!> For explicit methods, it directly uses the transport matrix X.
!>
!> Mathematical formulation for implicit methods:
!> \f[
!> \text{Solve: } A \cdot c^{n+1} = X \cdot c^n + Z \cdot c_{bd} + Y \cdot c_{rech} + f
!> \f]
!> where mixing ratios are extracted from \f$ A^{-1} \f$ and the right-hand side matrices.
!>
!> \param[in,out] this Transport object containing spatial/temporal discretization and BCs
!> \param[out] A_mat_lumped Optional lumped system matrix for mass balance computations [-]

subroutine compute_mixing_ratios_Delta_t_homog_2D(this)
    use BCs_subroutines_m
    use transport_transient_m, only: transport_2D_transient_c
    use arrays_m, only: diag_matrix_c,pentadiag_matrix_c, copy_real_array, inv_SPD_mat
    use time_discr_m, only: time_discr_homog_c
    use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c
    use penta_cholesky_inv_m, only: penta_spd_inverse
    implicit none
    
    class(transport_2D_transient_c) :: this
    
    integer(kind=4) :: i,j
    integer(kind=4) :: num_mix_rat
    integer(kind=4) :: num_inf_ext_mix_rat
    integer(kind=4) :: num_out_ext_mix_rat
    integer(kind=4) :: num_rech_mix_rat
    integer(kind=4) :: Nx, Ny
    integer(kind=4) :: row_i, col_i
    integer(kind=4) :: num_upstream, num_downstream
    integer(kind=4) :: nUpKept, nDownKept, dim_conc, pos  !> sparse-mixing bookkeeping (kept neighbours, shared by the concentration and reaction ratios)
    real(kind=8) :: theta,Delta_t
    real(kind=8), parameter :: tol_inv=1d-12
    real(kind=8), parameter :: tol_mix=1d-3  !> threshold: domain mixing ratios below this are dropped (sparse storage)
    real(kind=8) :: bd1, bd2, kept_dom_sum, dom_target_sum, scale  !> renormalisation of kept domain ratios
    real(kind=8), allocatable :: A_mat_inv(:,:),A_mat(:,:)
    real(kind=8) :: mix_sum  !> DEBUG: sum of mixing ratios for mass balance check
    real(kind=8) :: t1, t2  !> CPU timing variables
    
    type(pentadiag_matrix_c) :: E_mat,E_mat_prev,X_mat_T
    
!> Extract 2D grid dimensions
    select type (mesh=>this%spatial_discr)
    type is (mesh_2D_Euler_homog_c)
        Nx = mesh%Num_cells_x
        Ny = mesh%Num_cells_y
    class default
        error stop "2D homogeneous mesh expected in compute_mixing_ratios_Delta_t_homog_2D"
    end select
    
    write(*,*) ''
    write(*,*) '================================================================================'
    write(*,*) 'DEBUG [compute_mixing_ratios_Delta_t_homog_2D]: ENTRY'
    write(*,*) '================================================================================'
    write(*,*) '  Nx =', Nx, '  Ny =', Ny
    write(*,*) '  Num_targets =', this%spatial_discr%Num_targets
    write(*,*) '  Lagr_flag   =', this%Lagr_flag
    write(*,*) '  int_method  =', this%time_discr%int_method_tpt
    write(*,*) '  BCs labels  =', this%BCs%labels
    
!> Step 1: Allocate mixing ratio arrays for concentrations and reaction amounts
!> mixing_ratios_R shares the layout of mixing_ratios_conc / mix_conc_indices ([self, upstream...,
!> downstream...], externals included): the reaction ratio of a water that mixes its concentrations
!> but not its reaction rates (boundary and recharge waters) is stored as 0.
    call this%mixing_ratios_conc%allocate_array(this%spatial_discr%Num_targets)
    call this%mixing_ratios_R%allocate_array(this%spatial_discr%Num_targets)
    
!> Step 2: Determine number of boundary mixing ratios based on boundary condition types
    num_mix_rat=0
    num_inf_ext_mix_rat=1 !> Inflow from left boundary
    num_out_ext_mix_rat=1 !> Outflow to right boundary
    num_rech_mix_rat=0
    num_mix_rat=num_mix_rat+num_inf_ext_mix_rat+num_out_ext_mix_rat+num_rech_mix_rat
    
    write(*,*) '  num_inf_ext_mix_rat =', num_inf_ext_mix_rat
    write(*,*) '  num_out_ext_mix_rat =', num_out_ext_mix_rat
    write(*,*) '  num_rech_mix_rat    =', num_rech_mix_rat
    write(*,*) '  num_mix_rat (ext)   =', num_mix_rat
    
!> Step 3: Extract time step size and verify time discretization is homogeneous
    select type (time_discr=>this%time_discr)
    class is (time_discr_homog_c)
        Delta_t=time_discr%Delta_t
    class default
        error stop "This subroutine is only applied if time discretisation is homogeneous"
    end select

    write(*,*) '  Delta_t =', Delta_t

!> Step 4: Compute PDE matrices (transport, recharge, source, and mass matrices)
    call this%compute_trans_mat_PDE()
    call this%compute_rech_mat_PDE()
    call this%compute_source_term_PDE()
    call this%compute_F_mat_PDE()
    !> We impose the boundary conditions
    if (this%BCs%labels(1).eq.1 .and. this%BCs%labels(2).eq.1) then
        call Dirichlet_BCs_PDE(this)
    else if (this%BCs%labels(1).eq.1 .and. this%BCs%labels(2).eq.2) then
        call Dirichlet_Neumann_BCs_PDE(this)
    else if (this%BCs%labels(1).eq.2 .and. this%BCs%labels(2).eq.2) then
        call Neumann_homog_BCs(this)
    else if (this%BCs%labels(1).eq.3 .and. this%BCs%labels(2).eq.2) then
        call Robin_Neumann_homog_BCs(this)
    else if (this%BCs%labels(1).eq.3 .and. this%BCs%labels(2).eq.1) then
        call Robin_Dirichlet_BCs(this)
    else
        error stop "Boundary conditions not implemented yet"
    end if
!> Verify that rows of trans_mat sum to 0 (conservation check)
    do i = 1, this%spatial_discr%Num_targets
        col_i = mod(i-1, Nx) + 1
        row_i = (i-1)/Nx + 1
        mix_sum = this%trans_mat%diag(i)
        !> left neighbor (sub-diagonal, offset -1)
        if (col_i > 1) mix_sum = mix_sum + this%trans_mat%sub((row_i-1)*(Nx-1) + col_i - 1)
        !> right neighbor (super-diagonal, offset +1)
        if (col_i < Nx) mix_sum = mix_sum + this%trans_mat%super((row_i-1)*(Nx-1) + col_i)
        !> below neighbor (sub2-diagonal, offset -Nx)
        if (row_i > 1) mix_sum = mix_sum + this%trans_mat%sub2((col_i-1)*(Ny-1) + row_i - 1)
        !> above neighbor (super2-diagonal, offset +Nx)
        if (row_i < Ny) mix_sum = mix_sum + this%trans_mat%super2((col_i-1)*(Ny-1) + row_i)
        if (abs(mix_sum) > 1d-6) then
            write(*,'(A,I4,A,ES14.6)') 'ERROR: trans_mat row ', i, ' does not sum to 0: ', mix_sum
            error stop "trans_mat row does not sum to 0"
        end if
    end do
    if (this%time_discr%int_method_tpt.eq.1) then
        theta=0d0
!> Allocate dims based on actual 2D neighbors (left-to-right, bottom-to-top ordering)
        do i=1,this%mixing_ratios_conc%num_cols
            col_i = mod(i-1, Nx) + 1
            row_i = (i-1)/Nx + 1
            num_upstream = 0
            num_downstream = 0
            if (col_i > 1) num_upstream = num_upstream + 1    !> left neighbor
            if (row_i > 1) num_upstream = num_upstream + 1    !> below neighbor
            if (col_i < Nx) num_downstream = num_downstream + 1  !> right neighbor
            if (row_i < Ny) num_downstream = num_downstream + 1  !> above neighbor
            this%mixing_ratios_conc%cols(i)%dim = 1 + num_upstream + num_downstream + num_rech_mix_rat
            this%mixing_ratios_R%cols(i)%dim = this%mixing_ratios_conc%cols(i)%dim !> reaction ratios share the concentration layout
        end do
        call this%mixing_ratios_conc%allocate_columns()
        call this%mixing_ratios_R%allocate_columns()
    else
        num_mix_rat=num_mix_rat+this%spatial_discr%Num_targets
!> Memory budget check: implicit method requires O(N^2) storage for dense mixing ratios
        if (int(this%spatial_discr%Num_targets,8) * int(num_mix_rat,8) * 8_8 > 2000000000_8) then
            write(*,*) 'ERROR: Implicit method (int_method=', this%time_discr%int_method_tpt, &
                ') requires too much memory for this grid size.'
            write(*,*) '  Num_targets =', this%spatial_discr%Num_targets
            write(*,*) '  num_mix_rat =', num_mix_rat
            write(*,*) '  Estimated memory for mixing ratios: ', &
                int(this%spatial_discr%Num_targets,8) * int(num_mix_rat,8) * 8_8 / (1024_8*1024_8), ' MB'
            write(*,*) '  Consider using explicit method (int_method=1) for large 2D grids.'
            error stop "Implicit mixing ratios storage exceeds 2 GB memory budget"
        end if
        do i=1,this%mixing_ratios_conc%num_cols
            call this%mixing_ratios_conc%cols(i)%set_dim(num_mix_rat)
            call this%mixing_ratios_conc%cols(i)%allocate_vector()
            call this%mixing_ratios_R%cols(i)%set_dim(num_mix_rat) !> reaction ratios share the concentration layout (resized to the sparse kept-set below)
            call this%mixing_ratios_R%cols(i)%allocate_vector()
        end do
        if (this%time_discr%int_method_tpt.eq.2) then
            theta=1d0
        else if (this%time_discr%int_method_tpt.eq.3) then
            theta=5d-1
        else
            error stop "Time discretisation not implemented yet"
        end if
    end if
    call this%allocate_mix_conc_indices()

    write(*,*) '  theta              =', theta
    write(*,*) '  num_mix_rat (total)=', num_mix_rat
    write(*,*) '  mixing_ratios_conc%num_cols =', this%mixing_ratios_conc%num_cols
    write(*,*) '  mixing_ratios_R%num_cols    =', this%mixing_ratios_R%num_cols
    ! do i=1,this%mixing_ratios_conc%num_cols
    !     write(*,'(A,I4,A,I4,A,I4,A,I4)') '    target ', i, &
    !         ': conc_dim=', this%mixing_ratios_conc%cols(i)%dim, &
    !         '  R_dim=', this%mixing_ratios_R%cols(i)%dim, &
    !         '  idx_dim=', this%mix_conc_indices%cols(i)%dim
    ! end do

!> Step 6: Compute linear system matrices (E, X, A, Y, Z) and source vector (f)
    call this%compute_E_mat_2D(Delta_t,E_mat,E_mat_prev)
    call this%compute_X_mat_2D(theta,E_mat_prev)
    call this%compute_A_mat_2D(theta,E_mat)
    call this%compute_Y_mat()
    call this%compute_Z_mat()
    call this%compute_f_vec(this%time_discr%get_Delta_t())
    
    write(*,*) ''
    write(*,*) '--- PDE Matrices ---'
    ! write(*,*) '  F_mat%diag =', this%F_mat%diag
    ! write(*,*) '  trans_mat%diag  =', this%trans_mat%diag
    ! write(*,*) '  trans_mat%sub   =', this%trans_mat%sub
    ! write(*,*) '  trans_mat%super =', this%trans_mat%super
    ! write(*,*) '  trans_mat%sub2  =', this%trans_mat%sub2
    ! write(*,*) '  trans_mat%super2=', this%trans_mat%super2
    ! write(*,*) '  X_mat%diag  =', this%X_mat%diag
    ! write(*,*) '  X_mat%sub   =', this%X_mat%sub
    ! write(*,*) '  X_mat%super =', this%X_mat%super
    ! write(*,*) '  X_mat%sub2  =', this%X_mat%sub2
    ! write(*,*) '  X_mat%super2=', this%X_mat%super2
    ! write(*,*) '  A_mat%diag  =', this%A_mat%diag
    ! write(*,*) '  A_mat%sub   =', this%A_mat%sub
    ! write(*,*) '  A_mat%super =', this%A_mat%super
    ! write(*,*) '  A_mat%sub2  =', this%A_mat%sub2
    ! write(*,*) '  A_mat%super2=', this%A_mat%super2
    ! write(*,*) '  Y_mat%diag =', this%Y_mat%diag
    ! write(*,*) '  Z_mat      =', this%Z_mat
    ! write(*,*) '  f_vec      =', this%f_vec
    
!> Step 7: Compute mixing ratios from system matrices
    if (theta>0d0) then
!> Step 7a: Allocate arrays for implicit method and invert system matrix A
        call this%allocate_A_mat_inv()
        call this%allocate_mixing_ratios_mat_conc_mesh()
        call this%allocate_mixing_ratios_mat_conc_bd()
        if (this%Lagr_flag .and. all(this%BCs%labels>1)) then
            !A_mat=this%A_mat%fill_pentadiag_mat(this%spatial_discr%get_num_cells(1))
            !call inv_SPD_mat(A_mat,this%A_mat_inv)
            call cpu_time(t1)
            call penta_spd_inverse(this%spatial_discr%Num_targets, Nx, Ny, &
                this%A_mat%diag, this%A_mat%sub, this%A_mat%sub2, this%A_mat_inv)
            call cpu_time(t2)
            print *, 'CPU time of inverse of A_mat: ', t2 - t1, ' seconds'
        else
            write(*,*) ''
            write(*,*) '  WARNING: A_mat inversion path NOT taken (Lagr_flag=', this%Lagr_flag, ', BCs=', this%BCs%labels, ')'
            write(*,*) '  A_mat_inv is NOT computed — mixing ratios will be zero/uninitialized!'
        end if
        if (theta.eq.1d0) then
!> Step 7b: Fully implicit Euler - mixing ratios from A^{-1} only
            this%mixing_ratios_mat_conc_mesh=transpose(this%A_mat_inv)
                this%mixing_ratios_mat_conc_bd(:,1)=this%Z_mat(1)*this%mixing_ratios_mat_conc_mesh(1,:)
                this%mixing_ratios_mat_conc_bd(:,2)=this%Z_mat(2)*this%mixing_ratios_mat_conc_mesh(this%spatial_discr%Num_targets,:)
        else
!> Step 7c: Partially implicit - mixing ratios from X^T · A^{-1}
            call this%X_mat%compute_transpose_tridiag_matrix(X_mat_T)
            this%mixing_ratios_mat_conc_mesh=X_mat_T%prod_tridiag_mat_mat(transpose(this%A_mat_inv))
                this%mixing_ratios_mat_conc_bd(:,1)=&
                    (theta*this%Z_mat(1)+(1d0-theta)*this%Z_mat_prev(1))*this%A_mat_inv(:,1)
                this%mixing_ratios_mat_conc_bd(:,2)=&
                    (theta*this%Z_mat(2)+(1d0-theta)*this%Z_mat_prev(2))*this%A_mat_inv(:,this%spatial_discr%Num_targets)
        end if

        write(*,*) ''
        write(*,*) '--- mixing_ratios_mat_conc_mesh (domain mixing matrix) ---'
        ! do i=1,size(this%mixing_ratios_mat_conc_mesh,1)
        !     write(*,'(A,I4,A,*(ES14.6))') '    row ', i, ':', this%mixing_ratios_mat_conc_mesh(i,:)
        ! end do
        write(*,*) ''
        write(*,*) '--- mixing_ratios_mat_conc_bd (boundary mixing) ---'
        ! do i=1,size(this%mixing_ratios_mat_conc_bd,1)
        !     write(*,'(A,I4,A,*(ES14.6))') '    row ', i, ':', this%mixing_ratios_mat_conc_bd(i,:)
        ! end do

    !> Step 7e: Populate mixing ratios for each target (implicit methods) -- SPARSE storage.
    !> Only domain mixing ratios with |lambda| >= tol_mix are stored; the self ratio and the
    !> inflow/outflow boundary ratios are always kept. The kept domain ratios (self + retained
    !> neighbours) are renormalised so that, together with the exact boundary ratios, each column
    !> still sums to 1 (conservative mixing). This collapses the dense O(N^2) implicit mixing
    !> operator into a sparse one, drastically reducing setup and per-time-step mixing cost.
    !> The reaction-rate mixing ratios (mixing_ratios_R) are stored at the SAME positions as the
    !> concentration ones (mix_conc_indices layout, kept-set shared): A_mat_inv(i,j) for the kept
    !> domain neighbours, 0 at the boundary/recharge positions (they do not mix reaction rates).
        write(*,*) ''
        write(*,*) '--- Populating mix_conc_indices and mixing_ratios (implicit, sparse, tol=', tol_mix, ') ---'
        do i=1,this%spatial_discr%num_targets
            !> Column-major upstream/downstream classification (flow is along +x).
            !> A domain neighbour j is UPSTREAM iff its physical column (x) is strictly to the
            !> left of cell i's column; otherwise it is DOWNSTREAM. Same-column (same-x)
            !> cross-stream neighbours are therefore treated as downstream (NOT upstream),
            !> so only genuine x-upstream cells enter the implicit reaction term. The physical
            !> column of a (row-major) index is col = mod(index-1, Nx) + 1.
            col_i = mod(i-1, Nx) + 1
            !> Pass 1: count kept upstream (col<col_i) and downstream (col>=col_i) domain neighbours
            !> (|conc ratio| >= tol_mix).
            nUpKept=0
            nDownKept=0
            do j=1,this%spatial_discr%num_targets
                if (j==i) cycle
                if (mod(j-1, Nx) + 1 < col_i) then
                    !> upstream (physical column strictly to the left)
                    if (abs(this%mixing_ratios_mat_conc_mesh(j,i))>=tol_mix) nUpKept=nUpKept+1
                else
                    !> downstream (right column or same-x cross-stream)
                    if (abs(this%mixing_ratios_mat_conc_mesh(j,i))>=tol_mix) nDownKept=nDownKept+1
                end if
            end do
            !> Sparse column length: self + inflow ext + kept up + kept down + outflow ext + recharge
            dim_conc=1+num_inf_ext_mix_rat+nUpKept+nDownKept+num_out_ext_mix_rat+num_rech_mix_rat
            call this%mixing_ratios_conc%cols(i)%allocate_vector(dim_conc)
            call this%mix_conc_indices%cols(i)%allocate_vector(dim_conc+2)
            call this%mixing_ratios_R%cols(i)%allocate_vector(dim_conc)
            this%mixing_ratios_R%cols(i)%col_1=0d0 !> externals (boundary, recharge) do not mix reaction rates

            !> Self + boundary (kept exact) entries and the index/metadata bookkeeping
            bd1=this%mixing_ratios_mat_conc_bd(i,1)
            bd2=this%mixing_ratios_mat_conc_bd(i,2)
            this%mixing_ratios_conc%cols(i)%col_1(1)=this%mixing_ratios_mat_conc_mesh(i,i)             !> self (renormalised below)
            this%mixing_ratios_conc%cols(i)%col_1(1+num_inf_ext_mix_rat)=bd1                           !> inflow ext (exact)
            this%mixing_ratios_conc%cols(i)%col_1(dim_conc-num_rech_mix_rat)=bd2                        !> outflow ext (exact)
            this%mixing_ratios_R%cols(i)%col_1(1)=this%A_mat_inv(i,i)
            this%mix_conc_indices%cols(i)%col_1(1)=i+num_inf_ext_mix_rat
            this%mix_conc_indices%cols(i)%col_1(1+num_inf_ext_mix_rat)=num_inf_ext_mix_rat
            this%mix_conc_indices%cols(i)%col_1(dim_conc-num_rech_mix_rat)=&
                this%spatial_discr%num_targets+num_inf_ext_mix_rat+num_out_ext_mix_rat
            this%mix_conc_indices%cols(i)%col_1(this%mix_conc_indices%cols(i)%dim-1)=num_inf_ext_mix_rat+nUpKept
            this%mix_conc_indices%cols(i)%col_1(this%mix_conc_indices%cols(i)%dim)=nDownKept+num_out_ext_mix_rat

            !> Pass 2: fill kept domain neighbours (upstream then downstream) and accumulate kept sum.
            !> The reaction ratio of each kept neighbour goes to the same position "pos".
            kept_dom_sum=this%mixing_ratios_mat_conc_mesh(i,i)
            pos=1+num_inf_ext_mix_rat
!> Kept upstream targets (physical column < col_i)
            do j=1,this%spatial_discr%num_targets
                if (j==i) cycle
                if (mod(j-1, Nx) + 1 < col_i) then
                    if (abs(this%mixing_ratios_mat_conc_mesh(j,i))>=tol_mix) then
                        pos=pos+1
                        this%mixing_ratios_conc%cols(i)%col_1(pos)=this%mixing_ratios_mat_conc_mesh(j,i)
                        this%mixing_ratios_R%cols(i)%col_1(pos)=this%A_mat_inv(i,j)
                        this%mix_conc_indices%cols(i)%col_1(pos)=num_inf_ext_mix_rat+j
                        kept_dom_sum=kept_dom_sum+this%mixing_ratios_mat_conc_mesh(j,i)
                    end if
                end if
            end do
!> Kept downstream targets (physical column >= col_i, i.e. right column or same-x cross-stream)
            do j=1,this%spatial_discr%num_targets
                if (j==i) cycle
                if (mod(j-1, Nx) + 1 >= col_i) then
                    if (abs(this%mixing_ratios_mat_conc_mesh(j,i))>=tol_mix) then
                        pos=pos+1
                        this%mixing_ratios_conc%cols(i)%col_1(pos)=this%mixing_ratios_mat_conc_mesh(j,i)
                        this%mixing_ratios_R%cols(i)%col_1(pos)=this%A_mat_inv(i,j)
                        this%mix_conc_indices%cols(i)%col_1(pos)=num_inf_ext_mix_rat+j
                        kept_dom_sum=kept_dom_sum+this%mixing_ratios_mat_conc_mesh(j,i)
                    end if
                end if
            end do
!> Assign recharge mixing ratios (if any), appended after the outflow boundary entry
            do j=1,num_rech_mix_rat
                this%mixing_ratios_conc%cols(i)%col_1(dim_conc-num_rech_mix_rat+j)=&
                    this%Y_mat%diag(j)*this%A_mat_inv(i,j)
                this%mix_conc_indices%cols(i)%col_1(dim_conc-num_rech_mix_rat+j)=&
                    this%spatial_discr%num_targets+num_inf_ext_mix_rat+num_out_ext_mix_rat+j
            end do

!> Renormalise domain entries (self + kept neighbours) so the column sums to 1 with the
!> boundary (and recharge) ratios held exact. Discarded weight is redistributed proportionally.
            dom_target_sum=1d0-bd1-bd2
            do j=1,num_rech_mix_rat
                dom_target_sum=dom_target_sum-this%mixing_ratios_conc%cols(i)%col_1(dim_conc-num_rech_mix_rat+j)
            end do
            if (kept_dom_sum>0d0) then
                scale=dom_target_sum/kept_dom_sum
            else
                scale=1d0
            end if
            this%mixing_ratios_conc%cols(i)%col_1(1)=scale*this%mixing_ratios_conc%cols(i)%col_1(1)
            pos=1+num_inf_ext_mix_rat
            do j=1,nUpKept+nDownKept
                pos=pos+1
                this%mixing_ratios_conc%cols(i)%col_1(pos)=scale*this%mixing_ratios_conc%cols(i)%col_1(pos)
            end do

            !> Verify non-negativity and mass balance
            do j=1,this%mixing_ratios_conc%cols(i)%dim
                if (this%mixing_ratios_conc%cols(i)%col_1(j) < 0d0) then
                    write(*,'(A,I4,A,I4,A,ES15.7)') &
                        'WARNING: negative mixing ratio at target ', i, &
                        ', index ', j, ', value = ', this%mixing_ratios_conc%cols(i)%col_1(j)
                end if
            end do
            mix_sum = sum(this%mixing_ratios_conc%cols(i)%col_1)
            if (abs(mix_sum - 1d0) > 1d-5) then
                write(*,'(A,I4,A,ES14.6)') 'WARNING: mixing ratios column ', i, ' does not sum to 1: ', mix_sum
            end if
            ! write(*,'(A,I4,A)') '  --- Target i=', i, ' ---'
            ! write(*,'(A,*(ES12.4))') '    mixing_ratios_conc = ', this%mixing_ratios_conc%cols(i)%col_1
            ! write(*,'(A,ES14.6,A,ES10.2)') '    SUM(mixing_ratios_conc) = ', mix_sum, &
            !     '   error from 1 = ', abs(mix_sum - 1d0)
            ! write(*,'(A,*(ES12.4))') '    mixing_ratios_R    = ', this%mixing_ratios_R%cols(i)%col_1
            ! write(*,'(A,*(I6))') '    mix_conc_indices   = ', this%mix_conc_indices%cols(i)%col_1
            ! write(*,'(A,I4)') '      idx(1)   [current water]  = ', this%mix_conc_indices%cols(i)%col_1(1)
            ! write(*,'(A,I4)') '      idx(dim-1)[num upstream]  = ', &
            !     this%mix_conc_indices%cols(i)%col_1(this%mix_conc_indices%cols(i)%dim-1)
            ! write(*,'(A,I4)') '      idx(dim)  [num downstream]= ', &
            !     this%mix_conc_indices%cols(i)%col_1(this%mix_conc_indices%cols(i)%dim)
        end do
    else
!> Step 8: Assign mixing ratios directly from X matrix for explicit methods (2D)
!> Targets are ordered left-to-right, bottom-to-top in the 2D mesh
        write(*,*) ''
        write(*,*) '--- Populating mix_conc_indices and mixing_ratios (EXPLICIT 2D) ---'
        do i=1,this%spatial_discr%num_targets
            col_i = mod(i-1, Nx) + 1
            row_i = (i-1)/Nx + 1
            
            !> Self mixing ratio (diagonal)
            this%mixing_ratios_conc%cols(i)%col_1(1) = this%X_mat%diag(i)
            this%mix_conc_indices%cols(i)%col_1(1) = i + num_inf_ext_mix_rat
            this%mixing_ratios_R%cols(i)%col_1 = 0d0 !> explicit: reactions mix only with the target itself
            this%mixing_ratios_R%cols(i)%col_1(1) = 1d0
            
            num_upstream = 0
            num_downstream = 0
            j = 1  !> position counter (1 = self)
            
            !> Below neighbor (upstream, index i-Nx)
            if (row_i > 1) then
                j = j + 1
                num_upstream = num_upstream + 1
                this%mixing_ratios_conc%cols(i)%col_1(j) = this%X_mat%sub2((col_i-1)*(Ny-1) + row_i - 1)
                this%mix_conc_indices%cols(i)%col_1(j) = (i - Nx) + num_inf_ext_mix_rat
            end if
            !> Left neighbor (upstream, index i-1)
            if (col_i > 1) then
                j = j + 1
                num_upstream = num_upstream + 1
                this%mixing_ratios_conc%cols(i)%col_1(j) = this%X_mat%sub((row_i-1)*(Nx-1) + col_i - 1)
                this%mix_conc_indices%cols(i)%col_1(j) = (i - 1) + num_inf_ext_mix_rat
            end if
            !> Right neighbor (downstream, index i+1)
            if (col_i < Nx) then
                j = j + 1
                num_downstream = num_downstream + 1
                this%mixing_ratios_conc%cols(i)%col_1(j) = this%X_mat%super((row_i-1)*(Nx-1) + col_i)
                this%mix_conc_indices%cols(i)%col_1(j) = (i + 1) + num_inf_ext_mix_rat
            end if
            !> Above neighbor (downstream, index i+Nx)
            if (row_i < Ny) then
                j = j + 1
                num_downstream = num_downstream + 1
                this%mixing_ratios_conc%cols(i)%col_1(j) = this%X_mat%super2((col_i-1)*(Ny-1) + row_i)
                this%mix_conc_indices%cols(i)%col_1(j) = (i + Nx) + num_inf_ext_mix_rat
            end if
            
            !> Set upstream and downstream counts
            this%mix_conc_indices%cols(i)%col_1(this%mix_conc_indices%cols(i)%dim-1) = num_upstream
            this%mix_conc_indices%cols(i)%col_1(this%mix_conc_indices%cols(i)%dim) = num_downstream

            do j=1,this%mixing_ratios_conc%cols(i)%dim
                if (this%mixing_ratios_conc%cols(i)%col_1(j) < 0d0) then
                    write(*,'(A,I4,A,I4,A,ES15.7)') &
                        'WARNING: negative mixing ratio at target ', i, &
                        ', index ', j, ', value = ', this%mixing_ratios_conc%cols(i)%col_1(j)
                end if
            end do
            mix_sum = sum(this%mixing_ratios_conc%cols(i)%col_1)
            if (abs(mix_sum - 1d0) > 1d-5) then
                write(*,'(A,I4,A,ES14.6)') 'WARNING: mixing ratios column ', i, ' does not sum to 1: ', mix_sum
            end if
            ! write(*,'(A,I4,A,I4,A,I4,A)') '  --- Target i=', i, ' (col=', col_i, ', row=', row_i, ', explicit 2D) ---'
            ! write(*,'(A,*(ES12.4))') '    mixing_ratios_conc = ', this%mixing_ratios_conc%cols(i)%col_1
            ! write(*,'(A,ES14.6,A,ES10.2)') '    SUM = ', mix_sum, '   err = ', abs(mix_sum - 1d0)
            ! write(*,'(A,*(I6))') '    mix_conc_indices   = ', this%mix_conc_indices%cols(i)%col_1
        end do
    end if !> End of time integration method selection (implicit vs explicit)
    call copy_real_array(this%mixing_ratios_R, this%mixing_ratios_R_init)

    write(*,*) ''
    write(*,*) '--- Final mixing_ratios_R_init ---'
    ! do i=1,this%mixing_ratios_R_init%num_cols
    !     write(*,'(A,I4,A,*(ES12.4))') '  R_init(', i, ') = ', this%mixing_ratios_R_init%cols(i)%col_1
    ! end do

    write(*,*) ''
    write(*,*) '================================================================================'
    write(*,*) 'DEBUG [compute_mixing_ratios_Delta_t_homog_2D]: EXIT'
    write(*,*) '================================================================================'
    write(*,*) ''
end subroutine !> End of compute_mixing_ratios_Delta_t_homog subroutine 