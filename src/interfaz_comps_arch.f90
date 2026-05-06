!> This subroutine is the interface to solve a reactive mixing iteration for components
!> It has to be called only once in each time step
!> It is supposed to be called after the conservative transport has been solved
!> It uses Euler explicit and applies lumping to the kinetic mixing ratios
!> It reads the aqueous component concentrations after conservative transport in a file, and writes the concentrations after reactive mixing in a different file
subroutine interfaz_comps_arch(this,path,num_aq_comps,file_in,Delta_t,file_out)
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
    integer(kind=4) :: n_nc !> number of variable activity species
    integer(kind=4) :: n_v !> number of variable activity species per (single) reactive zone (for cv_hat)
    integer(kind=4) :: n_aq_eq,n_min_eq,n_exch_eq,n_gas_eq !> per-water equilibrium reaction counts
    integer(kind=4) :: n_aq_kin,n_min_kin !> per-water kinetic reaction counts
    integer(kind=4) :: n_eq_tot,n_kin_tot !> total eq/kin rates written per water
    logical :: CV_flag !> TRUE if converges, FALSE otherwise
    real(kind=8), allocatable :: u_tilde(:,:) !> concentration components after solving conservative transport
    real(kind=8), allocatable :: u_react(:) !> concentration components after solving reactive mixing
    real(kind=8), allocatable :: u_new(:,:) !> concentration components after solving reactive mixing
    real(kind=8), allocatable :: conc_nc(:,:) !> concentrations of variable activity species after solving reactive mixing
    real(kind=8), allocatable :: conc_nc_hat(:) !> conc var. act. species after transport mixing only (no reaction)
    real(kind=8), allocatable :: cv_hat(:) !> conc var. act. species after mixing+kinetic addition (input to compute_Re)
    real(kind=8), allocatable :: rk_used(:) !> kinetic reaction rates effectively used in this iteration (chem ordering)
    real(kind=8), allocatable :: Sk_nc(:,:) !> kinetic stoichiometric matrix (variable activity species)
    real(kind=8), allocatable :: Re_out(:,:) !> equilibrium reaction rates per target water [n_eq_tot, num_target_waters]
    real(kind=8), allocatable :: Rk_out(:,:) !> kinetic reaction rates per target water    [n_kin_tot, num_target_waters]
!> Pre-process
    allocate(u_tilde(num_aq_comps,this%num_target_waters),u_new(num_aq_comps,this%num_target_waters))
    allocate(u_react(num_aq_comps))
    n_nc=this%get_num_aq_var_act_species() !> we get number of variable activity species in the domain
    allocate(conc_nc(n_nc,this%num_target_waters)) !> we allocate conc_nc for all target waters in the domain
    !> Determine reaction-rate sizes from the first target water (assumed uniform across targets).
    associate(tw0 => this%waters(this%tar_wat_indices(1)))
        n_aq_eq  = tw0%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts
        n_min_eq = tw0%solid_chemistry%reactive_zone%num_minerals
        n_exch_eq= tw0%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
        if (associated(tw0%gas_chemistry)) then
            n_gas_eq = tw0%solid_chemistry%reactive_zone%gas_phase%num_gases_eq
        else
            n_gas_eq = 0
        end if
        n_aq_kin = tw0%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts + &
                   tw0%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        n_min_kin= tw0%solid_chemistry%mineral_zone%num_minerals_kin
    end associate
    n_eq_tot  = n_aq_eq + n_min_eq + n_exch_eq + n_gas_eq
    n_kin_tot = n_aq_kin + n_min_kin
    allocate(Re_out(max(n_eq_tot,1),this%num_target_waters))
    allocate(Rk_out(max(n_kin_tot,1),this%num_target_waters))
    Re_out=0d0
    Rk_out=0d0
!> Process
    !> We read the component concentrations after solving conservative transport
    open(unit=1,file=path//file_in,status='old',action='read')
    do i=1,num_aq_comps
        read(1,*) (u_tilde(i,j), j=1,this%num_target_waters) !> we read one row (component) at a time
    end do
    close(1)
    !> We solve reactive mixing for each target water
    do j=1,this%num_target_waters !> loop over target waters
        tw_idx = this%tar_wat_indices(j)
        associate(tw => this%waters(tw_idx))
        !> --------------------------------------------------------------------
        !> Step 1: pre-reaction speciation of u_tilde to obtain c_v_hat
        !> (variable activity species after transport mixing only). Needed as
        !> input to compute_Re for the equilibrium reaction rates below.
        !> --------------------------------------------------------------------
        if (n_eq_tot > 0) then
            allocate(conc_nc_hat(n_nc))
            call tw%compute_c_nc_from_u_Newton_ideal(tw%get_c1(),u_tilde(:,j),conc_nc_hat,niter,CV_flag)
            if (.not.CV_flag) then
                print *, "Target water index: ", tw_idx
                print *, "No convergence in pre-reaction speciation (u_tilde)"
                error stop
            end if
            cv_hat = tw%get_conc_nc() !> snapshot in canonical reactive-zone ordering
            !> Add kinetic contribution Δt·Sk_nc^T·rk_old to obtain the cv_hat
            !> expected by compute_Re (mixed + kinetic, before equilibrium step).
            if (n_kin_tot > 0) then
                rk_used = tw%get_rk_old()
                Sk_nc   = tw%get_Sk_nc()
                cv_hat  = cv_hat + Delta_t*matmul(transpose(Sk_nc),rk_used)
                deallocate(rk_used,Sk_nc)
            end if
            deallocate(conc_nc_hat)
        end if
        !> --------------------------------------------------------------------
        !> Step 2: chemical part of the components (kinetic contribution) and
        !> reactive mixing update of u.
        !> --------------------------------------------------------------------
        call tw%compute_react_term_EE_eq_kin(Delta_t,1.0d0,u_react) !> chemical part of components
        u_new(:,j)=u_tilde(:,j)+u_react !> we sum transport and reaction parts
        call tw%compute_c_nc_from_u_Newton_ideal(tw%get_c1(),u_new(:,j),conc_nc(:,j),niter,&
            CV_flag) !> we compute concentrations of variable activity species from component concentrations
        if (.not.CV_flag) then
            print *, "Target water index: ", tw_idx
            print *, "No convergence in speciation after reactive mixing iteration"
            print *, "Try reducing the time step"
            error stop
        end if
        !> --------------------------------------------------------------------
        !> Step 3: equilibrium reaction rates (compute_Re distributes Re into
        !> tw%Re, tw%solid_chemistry%Re and tw%gas_chemistry%Re) and copy to
        !> output buffer divided by Δt to express them as rates.
        !> --------------------------------------------------------------------
        if (n_eq_tot > 0) then
            call tw%compute_Re(cv_hat,Delta_t,1.0d0)
            deallocate(cv_hat)
            k=0
            do i=1,n_aq_eq
                Re_out(k+i,j) = tw%Re(i)/Delta_t
            end do
            k=k+n_aq_eq
            do i=1,n_min_eq
                Re_out(k+i,j) = tw%solid_chemistry%Re(i)/Delta_t
            end do
            k=k+n_min_eq
            do i=1,n_exch_eq
                Re_out(k+i,j) = tw%solid_chemistry%Re(n_min_eq+i)/Delta_t
            end do
            k=k+n_exch_eq
            if (n_gas_eq > 0) then
                do i=1,n_gas_eq
                    Re_out(k+i,j) = tw%gas_chemistry%Re(i)/Delta_t
                end do
            end if
        end if
        !> --------------------------------------------------------------------
        !> Step 4: kinetic reaction rates effectively used in this iteration
        !> (Euler-explicit uses rk_old). Order: aqueous (linear+redox) | mineral
        !> --------------------------------------------------------------------
        if (n_kin_tot > 0) then
            rk_used = tw%get_rk_old()
            do i=1,n_kin_tot
                Rk_out(i,j) = rk_used(i)
            end do
            deallocate(rk_used)
        end if
        end associate
    end do
!> Post-process
    !> We write the component concentrations after solving reactive mixing
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
    write(2,*) "Equilibrium reaction rates (rows: reactions [aqueous | minerals | exchange | gases], & 
        columns: targets):"
    write(2,*)
    if (n_eq_tot > 0) then
        do i=1,n_eq_tot
            write(2,"(2x,*(ES15.5))") (Re_out(i,j), j=1,this%num_target_waters)
        end do
    else
        write(2,*) " (no equilibrium reactions)"
    end if
    write(2,*)
    write(2,*) "Kinetic reaction rates (rows: reactions [aqueous (linear+redox) | minerals], & 
        columns: targets):"
    write(2,*)
    if (n_kin_tot > 0) then
        do i=1,n_kin_tot
            write(2,"(2x,*(ES15.5))") (Rk_out(i,j), j=1,this%num_target_waters)
        end do
    else
        write(2,*) " (no kinetic reactions)"
    end if
    close(2)
    deallocate(u_tilde,u_react,u_new,conc_nc,Re_out,Rk_out)
end subroutine
