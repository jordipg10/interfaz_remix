!> This subroutine is the interface to solve a reactive mixing iteration for components
!> when there are BOTH equilibrium and kinetic reactions.
!> It has to be called only once in each time step
!> It is supposed to be called after the conservative transport has been solved
!> It uses Euler explicit and applies lumping to the kinetic mixing ratios
!> It reads the aqueous component concentrations after conservative transport in a file, and writes the concentrations after reactive mixing in a different file
subroutine interfaz_comps_arch_eq_kin(this,path,num_aq_comps,file_in,Delta_t,file_out)
    use chemistry_m, only: chemistry_c
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: path !> path for input and output files
    integer(kind=4), intent(in) :: num_aq_comps !> number of aqueous components
    character(len=*), intent(in) :: file_in !> name of file containing aqueous component concentrations after solving conservative transport iteration
    real(kind=8), intent(in) :: Delta_t !> time step
    character(len=*), intent(in) :: file_out !> name of file containing variable activity species and aqueous component concentrations after solving reactive mixing iteration
!> Variables
    integer(kind=4) :: i,j,k,tw_idx !> loop variables / shorthand for target water index
    integer(kind=4) :: niter !> number of iterations in Newton algorithm
    integer(kind=4) :: n_nc !> number of variable activity species (full system: aq + exchange + gas var-act)
    integer(kind=4) :: n_p !> number of primary species per reactive zone
    integer(kind=4) :: n_aq_kin,n_min_kin,n_gas_kin !> per-water kinetic reaction counts
    integer(kind=4) :: n_kin_tot !> total kinetic rate rows written per water
    logical :: CV_flag !> TRUE if converges, FALSE otherwise
    real(kind=8), allocatable :: u_tilde(:,:) !> concentration components after solving conservative transport
    real(kind=8), allocatable :: u_react(:) !> chemical reaction contribution to component concentrations
    real(kind=8), allocatable :: u_new(:,:) !> concentration components after solving reactive mixing
    real(kind=8), allocatable :: conc_nc(:,:) !> concentrations of variable activity species after solving reactive mixing
    real(kind=8), allocatable :: rk_out(:,:) !> kinetic reaction rates per target water [n_kin_tot, num_target_waters]
!> Pre-process
    allocate(u_tilde(num_aq_comps,this%num_target_waters),u_new(num_aq_comps,this%num_target_waters))
    allocate(u_react(num_aq_comps))
    !> Determine reaction counts from the first target water (assumed uniform across targets).
    associate(tw0 => this%waters(this%tar_wat_indices(1)))
        n_aq_kin = tw0%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts + &
                   tw0%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        n_min_kin= tw0%solid_chemistry%mineral_zone%num_minerals_kin
        n_gas_kin= tw0%solid_chemistry%reactive_zone%chem_syst%num_gas_kin_reacts
        n_p      = tw0%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
        n_nc     = tw0%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    end associate
    allocate(conc_nc(n_nc,this%num_target_waters))
    n_kin_tot = n_aq_kin + n_min_kin + n_gas_kin
    allocate(rk_out(n_kin_tot,this%num_target_waters))
    rk_out=0d0
!> Process
    !> We read the component concentrations after solving conservative transport
    open(unit=1,file=path//file_in,status='old',action='read')
    do i=1,num_aq_comps
        read(1,*) (u_tilde(i,j), j=1,this%num_target_waters)
    end do
    close(1)
    !> We solve reactive mixing for each target water
    do j=1,this%num_target_waters
        tw_idx = this%tar_wat_indices(j)
        associate(tw => this%waters(tw_idx))
        !> --------------------------------------------------------------------
        !> Step 1: kinetic chemical contribution to components and reactive mixing
        !> update of u; then speciate u_new to obtain conc_nc.
        !> --------------------------------------------------------------------
        call tw%compute_react_term_EE_eq_kin(Delta_t,1.0d0,u_react)
        u_new(:,j)=u_tilde(:,j)+u_react
        u_new(:,j)=max(u_new(:,j),0d0) !> physical lower bound: total component concentrations cannot be negative (explicit-Euler overshoot clamp)
        call tw%compute_c_nc_from_u_Newton_ideal(tw%get_c1(),u_new(:,j),conc_nc(:,j),niter,CV_flag)
        if (.not.CV_flag) then
            print *, "Target water index: ", tw_idx
            print *, "No convergence in speciation after reactive mixing iteration"
            print *, "Try reducing the time step"
            error stop
        end if
        !> --------------------------------------------------------------------
        !> Step 2: kinetic reaction rates (Euler-explicit uses rk_old).
        !> Order: aqueous (linear+redox) | mineral | gas
        !> --------------------------------------------------------------------
        k=0
        do i=1,n_aq_kin
            rk_out(k+i,j) = tw%rk_old(i)
        end do
        k=k+n_aq_kin
        do i=1,n_min_kin
            rk_out(k+i,j) = tw%solid_chemistry%rk_old(i)
        end do
        !> Gas kinetic rates not yet implemented in gas_chemistry_c
        end associate
    end do
!> Post-process
    open(unit=2,file=path//file_out,status='unknown',form='formatted')
    write(2,*) "Aqueous component concentrations after solving reactive mixing iteration (rows: components, columns: & 
        targets):"
    write(2,*)
    do i=1,num_aq_comps
        write(2,"(2x,*(ES15.5))") (u_new(i,j), j=1,this%num_target_waters)
    end do
    write(2,*)
    write(2,*) "Concentrations of variable activity species after solving reactive mixing iteration (rows: species, & 
        columns: targets):"
    write(2,*)
    do i=1,n_nc
        write(2,"(2x,*(ES15.5))") (conc_nc(i,j), j=1,this%num_target_waters)
    end do
    write(2,*)
    write(2,*) "Kinetic reaction rates (rows: reactions [aqueous (linear+redox) | minerals], & 
        columns: targets):"
    write(2,*)
    do i=1,n_kin_tot
        write(2,"(2x,*(ES15.5))") (rk_out(i,j), j=1,this%num_target_waters)
    end do
    close(2)
    deallocate(u_tilde,u_react,u_new,conc_nc,rk_out)
end subroutine

!> This subroutine is the interface to solve a reactive mixing iteration for components
!> when there are ONLY equilibrium reactions (no kinetics).
!> It has to be called only once in each time step
!> It is supposed to be called after the conservative transport has been solved
!> Since there are no kinetic reactions, u_new == u_tilde and conc_nc is its speciation.
!> It reads the aqueous component concentrations after conservative transport in a file, and writes the concentrations after reactive mixing in a different file
subroutine interfaz_comps_arch_eq(this,path,num_aq_comps,file_in,Delta_t,file_out)
    use chemistry_m, only: chemistry_c
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: path !> path for input and output files
    integer(kind=4), intent(in) :: num_aq_comps !> number of aqueous components
    character(len=*), intent(in) :: file_in !> name of file containing aqueous component concentrations after solving conservative transport iteration
    real(kind=8), intent(in) :: Delta_t !> time step
    character(len=*), intent(in) :: file_out !> name of file containing variable activity species and aqueous component concentrations after solving reactive mixing iteration
!> Variables
    integer(kind=4) :: i,j,tw_idx !> loop variables / shorthand for target water index
    integer(kind=4) :: niter !> number of iterations in Newton algorithm
    integer(kind=4) :: n_nc !> number of variable activity species (full system: aq + exchange + gas var-act)
    logical :: CV_flag !> TRUE if converges, FALSE otherwise
    real(kind=8), allocatable :: u_tilde(:,:) !> concentration components after solving conservative transport (= u_new, no kinetics)
    real(kind=8), allocatable :: conc_nc(:,:) !> concentrations of variable activity species after solving reactive mixing
!> Pre-process
    allocate(u_tilde(num_aq_comps,this%num_target_waters))
    !> Determine n_nc from the first target water (assumed uniform across targets).
    associate(tw0 => this%waters(this%tar_wat_indices(1)))
        n_nc = tw0%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    end associate
    allocate(conc_nc(n_nc,this%num_target_waters))
!> Process
    !> We read the component concentrations after solving conservative transport
    open(unit=1,file=path//file_in,status='old',action='read')
    do i=1,num_aq_comps
        read(1,*) (u_tilde(i,j), j=1,this%num_target_waters)
    end do
    close(1)
    !> We solve reactive mixing for each target water
    !> (u_new == u_tilde since there are no kinetic reactions)
    do j=1,this%num_target_waters
        tw_idx = this%tar_wat_indices(j)
        associate(tw => this%waters(tw_idx))
        !> Speciate u_tilde to obtain conc_nc (new equilibrium state).
        call tw%compute_c_nc_from_u_Newton_ideal(tw%get_c1(),u_tilde(:,j),conc_nc(:,j),niter,CV_flag)
        if (.not.CV_flag) then
            print *, "Target water index: ", tw_idx
            print *, "No convergence in speciation after reactive mixing iteration"
            print *, "Try reducing the time step"
            error stop
        end if
        end associate
    end do
!> Post-process
    open(unit=2,file=path//file_out,status='unknown',form='formatted')
    write(2,*) "Aqueous component concentrations after solving reactive mixing iteration (rows: components, columns: & 
        targets):"
    write(2,*)
    do i=1,num_aq_comps
        write(2,"(2x,*(ES15.5))") (u_tilde(i,j), j=1,this%num_target_waters)
    end do
    write(2,*)
    write(2,*) "Concentrations of variable activity species after solving reactive mixing iteration (rows: species, & 
        columns: targets):"
    write(2,*)
    do i=1,n_nc
        write(2,"(2x,*(ES15.5))") (conc_nc(i,j), j=1,this%num_target_waters)
    end do
    write(2,*)
    close(2)
    deallocate(u_tilde,conc_nc)
end subroutine

!> This subroutine is the interface to solve a reactive mixing iteration for components
!> when there are BOTH equilibrium and kinetic reactions.
!> It has to be called only once in each time step
!> It is supposed to be called after the conservative transport has been solved
!> It uses Euler explicit and applies lumping to the kinetic mixing ratios
!> It reads the aqueous component concentrations after conservative transport in a file, and writes the concentrations after reactive mixing in a different file
!> In file_in, the rows are targets and columns are components (transposed w.r.t. interfaz_comps_arch_eq_kin)
subroutine interfaz_comps_arch_eq_kin_T(this,path,num_aq_comps,file_in,Delta_t,file_out)
    use chemistry_m, only: chemistry_c
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: path !> path for input and output files
    integer(kind=4), intent(in) :: num_aq_comps !> number of aqueous components
    character(len=*), intent(in) :: file_in !> name of file containing aqueous component concentrations after solving conservative transport iteration
    real(kind=8), intent(in) :: Delta_t !> time step
    character(len=*), intent(in) :: file_out !> name of file containing variable activity species and aqueous component concentrations after solving reactive mixing iteration
!> Variables
    integer(kind=4) :: i,j,k,tw_idx !> loop variables / shorthand for target water index
    integer(kind=4) :: niter !> number of iterations in Newton algorithm
    integer(kind=4) :: n_nc !> number of variable activity species (full system: aq + exchange + gas var-act)
    integer(kind=4) :: n_p !> number of primary species per reactive zone
    integer(kind=4) :: n_aq_kin,n_min_kin,n_gas_kin !> per-water kinetic reaction counts
    integer(kind=4) :: n_kin_tot !> total kinetic rate rows written per water
    logical :: CV_flag !> TRUE if converges, FALSE otherwise
    real(kind=8), allocatable :: u_tilde(:,:) !> concentration components after solving conservative transport
    real(kind=8), allocatable :: u_react(:) !> chemical reaction contribution to component concentrations
    real(kind=8), allocatable :: u_new(:,:) !> concentration components after solving reactive mixing
    real(kind=8), allocatable :: conc_nc(:,:) !> concentrations of variable activity species after solving reactive mixing
    real(kind=8), allocatable :: rk_out(:,:) !> kinetic reaction rates per target water [n_kin_tot, num_target_waters]
!> Pre-process
    allocate(u_tilde(num_aq_comps,this%num_target_waters),u_new(num_aq_comps,this%num_target_waters))
    allocate(u_react(num_aq_comps))
    !> Determine reaction counts from the first target water (assumed uniform across targets).
    associate(tw0 => this%waters(this%tar_wat_indices(1)))
        n_aq_kin = tw0%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts + &
                   tw0%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        n_min_kin= tw0%solid_chemistry%mineral_zone%num_minerals_kin
        n_gas_kin= tw0%solid_chemistry%reactive_zone%chem_syst%num_gas_kin_reacts
        n_p      = tw0%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
        n_nc     = tw0%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    end associate
    allocate(conc_nc(n_nc,this%num_target_waters))
    n_kin_tot = n_aq_kin + n_min_kin + n_gas_kin
    allocate(rk_out(n_kin_tot,this%num_target_waters))
    rk_out=0d0
!> Process
    !> We read the component concentrations after solving conservative transport
    !> (rows: targets, columns: components)
    open(unit=1,file=path//file_in,status='old',action='read')
    do j=1,this%num_target_waters
        read(1,*) (u_tilde(i,j), i=1,num_aq_comps)
    end do
    close(1)
    !> We solve reactive mixing for each target water
    do j=1,this%num_target_waters
        tw_idx = this%tar_wat_indices(j)
        associate(tw => this%waters(tw_idx))
        !> --------------------------------------------------------------------
        !> Step 1: kinetic chemical contribution to components and reactive mixing
        !> update of u; then speciate u_new to obtain conc_nc.
        !> --------------------------------------------------------------------
        call tw%compute_react_term_EE_eq_kin(Delta_t,1.0d0,u_react)
        u_new(:,j)=u_tilde(:,j)+u_react
        u_new(:,j)=max(u_new(:,j),0d0) !> physical lower bound: total component concentrations cannot be negative (explicit-Euler overshoot clamp)
        call tw%compute_c_nc_from_u_Newton_ideal(tw%get_c1(),u_new(:,j),conc_nc(:,j),niter,CV_flag)
        if (.not.CV_flag) then
            print *, "Target water index: ", tw_idx
            print *, "No convergence in speciation after reactive mixing iteration"
            print *, "Try reducing the time step"
            error stop
        end if
        !> --------------------------------------------------------------------
        !> Step 2: kinetic reaction rates (Euler-explicit uses rk_old).
        !> Order: aqueous (linear+redox) | mineral | gas
        !> --------------------------------------------------------------------
        k=0
        do i=1,n_aq_kin
            rk_out(k+i,j) = tw%rk_old(i)
        end do
        k=k+n_aq_kin
        do i=1,n_min_kin
            rk_out(k+i,j) = tw%solid_chemistry%rk_old(i)
        end do
        !> Gas kinetic rates not yet implemented in gas_chemistry_c
        end associate
    end do
!> Post-process
    open(unit=2,file=path//file_out,status='unknown',form='formatted')
    write(2,*) "Aqueous component concentrations after solving reactive mixing iteration (rows: components, columns: & 
        targets):"
    write(2,*)
    do i=1,num_aq_comps
        write(2,"(2x,*(ES15.5))") (u_new(i,j), j=1,this%num_target_waters)
    end do
    write(2,*)
    write(2,*) "Concentrations of variable activity species after solving reactive mixing iteration (rows: species, & 
        columns: targets):"
    write(2,*)
    do i=1,n_nc
        write(2,"(2x,*(ES15.5))") (conc_nc(i,j), j=1,this%num_target_waters)
    end do
    write(2,*)
    write(2,*) "Kinetic reaction rates (rows: reactions [aqueous (linear+redox) | minerals], & 
        columns: targets):"
    write(2,*)
    do i=1,n_kin_tot
        write(2,"(2x,*(ES15.5))") (rk_out(i,j), j=1,this%num_target_waters)
    end do
    close(2)
    deallocate(u_tilde,u_react,u_new,conc_nc,rk_out)
end subroutine

!> This subroutine is the interface to solve a reactive mixing iteration for components
!> when there are ONLY equilibrium reactions (no kinetics).
!> It has to be called only once in each time step
!> It is supposed to be called after the conservative transport has been solved
!> Since there are no kinetic reactions, u_new == u_tilde and conc_nc is its speciation.
!> It reads the aqueous component concentrations after conservative transport in a file, and writes the concentrations after reactive mixing in a different file
!> In file_in, the rows are targets and columns are components (transposed w.r.t. interfaz_comps_arch_eq)
subroutine interfaz_comps_arch_eq_T(this,path,num_aq_comps,file_in,Delta_t,file_out)
    use chemistry_m, only: chemistry_c
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: path !> path for input and output files
    integer(kind=4), intent(in) :: num_aq_comps !> number of aqueous components
    character(len=*), intent(in) :: file_in !> name of file containing aqueous component concentrations after solving conservative transport iteration
    real(kind=8), intent(in) :: Delta_t !> time step
    character(len=*), intent(in) :: file_out !> name of file containing variable activity species and aqueous component concentrations after solving reactive mixing iteration
!> Variables
    integer(kind=4) :: i,j,tw_idx !> loop variables / shorthand for target water index
    integer(kind=4) :: niter !> number of iterations in Newton algorithm
    integer(kind=4) :: n_nc !> number of variable activity species (full system: aq + exchange + gas var-act)
    logical :: CV_flag !> TRUE if converges, FALSE otherwise
    real(kind=8), allocatable :: u_tilde(:,:) !> concentration components after solving conservative transport (= u_new, no kinetics)
    real(kind=8), allocatable :: conc_nc(:,:) !> concentrations of variable activity species after solving reactive mixing
!> Pre-process
    allocate(u_tilde(num_aq_comps,this%num_target_waters))
    !> Determine n_nc from the first target water (assumed uniform across targets).
    associate(tw0 => this%waters(this%tar_wat_indices(1)))
        n_nc = tw0%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    end associate
    allocate(conc_nc(n_nc,this%num_target_waters))
!> Process
    !> We read the component concentrations after solving conservative transport
    !> (rows: targets, columns: components)
    open(unit=1,file=path//file_in,status='old',action='read')
    do j=1,this%num_target_waters
        read(1,*) (u_tilde(i,j), i=1,num_aq_comps)
    end do
    close(1)
    !> We solve reactive mixing for each target water
    !> (u_new == u_tilde since there are no kinetic reactions)
    do j=1,this%num_target_waters
        tw_idx = this%tar_wat_indices(j)
        associate(tw => this%waters(tw_idx))
        !> Speciate u_tilde to obtain conc_nc (new equilibrium state).
        call tw%compute_c_nc_from_u_Newton_ideal(tw%get_c1(),u_tilde(:,j),conc_nc(:,j),niter,CV_flag)
        if (.not.CV_flag) then
            print *, "Target water index: ", tw_idx
            print *, "No convergence in speciation after reactive mixing iteration"
            print *, "Try reducing the time step"
            error stop
        end if
        end associate
    end do
!> Post-process
    open(unit=2,file=path//file_out,status='unknown',form='formatted')
    write(2,*) "Aqueous component concentrations after solving reactive mixing iteration (rows: components, columns: & 
        targets):"
    write(2,*)
    do i=1,num_aq_comps
        write(2,"(2x,*(ES15.5))") (u_tilde(i,j), j=1,this%num_target_waters)
    end do
    write(2,*)
    write(2,*) "Concentrations of variable activity species after solving reactive mixing iteration (rows: species, & 
        columns: targets):"
    write(2,*)
    do i=1,n_nc
        write(2,"(2x,*(ES15.5))") (conc_nc(i,j), j=1,this%num_target_waters)
    end do
    write(2,*)
    close(2)
    deallocate(u_tilde,conc_nc)
end subroutine
