!> This subroutine is the interface to solve a FULL reactive mixing iteration for
!> components (BOTH the conservative mixing and the reactive mixing) when there are
!> equilibrium reactions, optionally combined with kinetic reactions.
!> It has to be called only once in each time step.
!> Unlike interfaz_comps_arch_eq_kin (which reads the component concentrations AFTER
!> the conservative transport has been solved externally), this interface reads the
!> component concentrations PRIOR to mixing together with the mixing ratios, performs
!> the conservative mixing internally, and only then solves the reactive mixing.
!> It uses Euler explicit and applies lumping to the kinetic mixing ratios.
!> Input files (both located in "path"):
!>   - file_in : aqueous component concentrations PRIOR to mixing.
!>               Rows: components, columns: waters (one column per mixing water).
!>   - file_mix: mixing ratios stored in row-major order. Rows: target waters,
!>               columns: waters. Entry (j,k) is the mixing ratio of water k into
!>               target water j, so that the conservatively mixed component
!>               concentration of target j is sum_k mix_ratios(j,k)*conc_prior(:,k).
!> The mixing matrix need NOT be square: the number of mixing waters (columns) may
!> differ from the number of target waters (rows). The number of mixing waters is
!> auto-detected from the number of columns in the first row of file_mix.
!> IMPORTANT: the column ordering (waters) of file_in must match that of file_mix,
!> i.e. file_in must have exactly as many columns as file_mix.
!> It writes the variable activity species and aqueous component concentrations after
!> solving the reactive mixing iteration in file_out.
subroutine interfaz_comps_arch_eq_kin_mix(this,path,num_aq_comps,file_in,file_mix,Delta_t,file_out)
    use chemistry_m, only: chemistry_c
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: path !> path for input and output files
    integer(kind=4), intent(in) :: num_aq_comps !> number of aqueous components
    character(len=*), intent(in) :: file_in !> name of file containing aqueous component concentrations PRIOR to mixing (rows: components, columns: waters)
    character(len=*), intent(in) :: file_mix !> name of file containing the mixing ratios in row-major order (rows: target waters, columns: waters)
    real(kind=8), intent(in) :: Delta_t !> time step
    character(len=*), intent(in) :: file_out !> name of file containing variable activity species and aqueous component concentrations after solving reactive mixing iteration
!> Variables
    integer(kind=4) :: i,j,k,tw_idx !> loop variables / shorthand for target water index
    integer(kind=4) :: niter !> number of iterations in Newton algorithm
    integer(kind=4) :: n_nc !> number of variable activity species (full system: aq + exchange + gas var-act)
    integer(kind=4) :: n_p !> number of primary species per reactive zone
    integer(kind=4) :: n_aq_kin,n_min_kin,n_gas_kin !> per-water kinetic reaction counts
    integer(kind=4) :: n_kin_tot !> total kinetic rate rows written per water
    integer(kind=4) :: num_mix_wat !> number of mixing waters (columns of the mixing matrix); may differ from num_target_waters
    integer(kind=4) :: ios_local !> I/O status for the column-count probe read
    character(len=:), allocatable :: first_line !> first row of file_mix, used to auto-detect num_mix_wat
    logical :: CV_flag !> TRUE if converges, FALSE otherwise
    real(kind=8), allocatable :: conc_prior(:,:) !> aqueous component concentrations PRIOR to mixing [num_aq_comps, num_mix_wat]
    real(kind=8), allocatable :: mix_ratios(:,:) !> mixing ratios [num_target_waters, num_mix_wat]; mix_ratios(j,k)=ratio of water k into target j
    real(kind=8), allocatable :: u_tilde(:,:) !> concentration components after conservative mixing
    real(kind=8), allocatable :: u_react(:) !> chemical reaction contribution to component concentrations
    real(kind=8), allocatable :: u_new(:,:) !> concentration components after solving reactive mixing
    real(kind=8), allocatable :: conc_nc(:,:) !> concentrations of variable activity species after solving reactive mixing
    real(kind=8), allocatable :: rk_out(:,:) !> kinetic reaction rates per target water [n_kin_tot, num_target_waters]
    character(len=256), allocatable :: react_names(:) !> names of kinetic reactions, in the same row order as rk_out
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
    !> Build reaction-name labels (aqueous linear, then aqueous redox, then
    !> kinetic minerals, then gas kinetic). Taken from the first target water
    !> since the chemical system is shared across all target waters.
    !> Only needed when the system actually has kinetic reactions.
    allocate(react_names(n_kin_tot))
    if (n_kin_tot > 0) then
        associate(tw0 => this%waters(this%tar_wat_indices(1)))
            k=0
            do i=1,tw0%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
                react_names(k+i)=tw0%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%name
            end do
            k=k+tw0%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
            do i=1,tw0%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
                react_names(k+i)=tw0%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%name
            end do
            k=k+tw0%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
            do i=1,n_min_kin
                react_names(k+i)=tw0%solid_chemistry%reactive_zone%chem_syst%minerals( &
                    tw0%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%mineral%name
            end do
            k=k+n_min_kin
            do i=1,n_gas_kin
                write(react_names(k+i),'("gas_kin_",i0)') i
            end do
        end associate
    end if
!> Process
    !> We read the mixing ratios in row-major order
    !> (rows: target waters, columns: waters). The number of mixing waters
    !> (columns) is auto-detected from the first row, so the matrix need not be
    !> square.
    open(unit=1,file=path//file_mix,status='old',action='read')
    call read_long_line(1,first_line,ios_local)
    if (ios_local /= 0) then
        print *, "Error reading the first row of the mixing-ratios file: ", path//file_mix
        error stop
    end if
    num_mix_wat = count_tokens(first_line)
    if (num_mix_wat < 1) then
        print *, "Could not detect any mixing-ratio columns in: ", path//file_mix
        error stop
    end if
    !> The rows of the mixing matrix are the target waters and its columns are the
    !> mixing waters. Every target water must itself be one of the mixing waters
    !> (a water always mixes with itself), so the number of mixing waters
    !> (columns) can never be smaller than the number of target waters (rows).
    !> Equivalently: the mixing matrix may have fewer rows than columns, but never
    !> more rows than columns.
    if (num_mix_wat < this%num_target_waters) then
        print *, "Invalid mixing-ratios matrix: more rows (target waters) than columns (mixing waters)."
        print *, "  Number of target waters (rows):    ", this%num_target_waters
        print *, "  Number of mixing waters (columns): ", num_mix_wat
        print *, "Each target water must be among the mixing waters, so columns >= rows is required."
        error stop
    end if
    rewind(1)
    allocate(mix_ratios(this%num_target_waters,num_mix_wat))
    do j=1,this%num_target_waters
        read(1,*) (mix_ratios(j,k), k=1,num_mix_wat)
    end do
    close(1)
    !> We read the aqueous component concentrations PRIOR to mixing
    !> (rows: components, columns: waters). The number of columns must match the
    !> number of mixing waters detected from the mixing-ratios file.
    allocate(conc_prior(num_aq_comps,num_mix_wat))
    open(unit=1,file=path//file_in,status='old',action='read')
    do i=1,num_aq_comps
        read(1,*) (conc_prior(i,j), j=1,num_mix_wat)
    end do
    close(1)
    !> --------------------------------------------------------------------
    !> Conservative mixing of the component concentrations:
    !>   u_tilde(:,j) = sum_k mix_ratios(j,k) * conc_prior(:,k)
    !> Components are linear combinations of species, so mixing them directly
    !> is equivalent to conservatively mixing the underlying species.
    !> Shapes: (num_aq_comps x num_mix_wat) * (num_mix_wat x num_target_waters)
    !>         = (num_aq_comps x num_target_waters).
    !> --------------------------------------------------------------------
    u_tilde = matmul(conc_prior, transpose(mix_ratios))
    !> We solve reactive mixing for each target water
    do j=1,this%num_target_waters
        tw_idx = this%tar_wat_indices(j)
        associate(tw => this%waters(tw_idx))
        if (n_kin_tot > 0) then
            !> ----------------------------------------------------------------
            !> Step 1: kinetic chemical contribution to components and reactive
            !> mixing update of u; then speciate u_new to obtain conc_nc.
            !> ----------------------------------------------------------------
            call tw%compute_react_term_EE_eq_kin(Delta_t,1.0d0,u_react)
            u_new(:,j)=u_tilde(:,j)+u_react
            !> Components are linear combinations of species and may legitimately
            !> be negative. However, a component whose row in comp_mat_aq is a
            !> single +1 entry IS a species concentration and must stay
            !> non-negative; clamp only those (explicit-Euler overshoot guard).
            block
                integer(kind=4) :: ic, nz_idx, nz_cnt
                real(kind=8)    :: cval
                do ic=1,num_aq_comps
                    nz_cnt=0; nz_idx=0
                    do k=1,size(tw%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq,2)
                        cval=tw%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq(ic,k)
                        if (cval/=0d0) then
                            nz_cnt=nz_cnt+1
                            nz_idx=k
                            if (nz_cnt>1) exit
                        end if
                    end do
                    if (nz_cnt==1 .and. &
                        tw%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq(ic,nz_idx)==1d0) then
                        if (u_new(ic,j)<0d0) u_new(ic,j)=0d0
                    end if
                end do
            end block
            call tw%compute_c_nc_from_u_Newton_ideal(tw%get_c1(),u_new(:,j),conc_nc(:,j),niter,CV_flag)
            if (.not.CV_flag) then
                print *, "Target water index: ", tw_idx
                print *, "No convergence in speciation after reactive mixing iteration"
                print *, "Try reducing the time step, or solve the reactive mixing implicitly."
                error stop
            end if
            !> Make u_new consistent with the speciated state: recompute it from
            !> the aqueous variable-activity species concentrations using
            !> comp_mat_aq. Without this step, a species-equivalent component that
            !> was clamped to 0 above can disagree with the small positive species
            !> concentration recovered by Newton (e.g. ~1e-15 for dissolved O2).
            u_new(:,j) = matmul( &
                tw%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq, &
                conc_nc(1:tw%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species, j))
            !> Sanity check: kinetic reactions must not drive any species
            !> concentration negative. If they do, the time step is too large.
            if (any(conc_nc(:,j) < 0d0)) then
                print *, "Target water index: ", tw_idx
                do i=1,n_nc
                    if (conc_nc(i,j) < 0d0) then
                        print *, "  Species index ", i, " concentration: ", conc_nc(i,j)
                    end if
                end do
                print *, "Kinetic reactions produced a negative species concentration."
                print *, "Try reducing the time step, or solve the reactive mixing implicitly."
                error stop
            end if
            !> ----------------------------------------------------------------
            !> Step 2: kinetic reaction rates (Euler-explicit uses rk_old).
            !> Order: aqueous (linear+redox) | mineral | gas
            !> ----------------------------------------------------------------
            k=0
            do i=1,n_aq_kin
                rk_out(k+i,j) = tw%rk_old(i)
            end do
            k=k+n_aq_kin
            do i=1,n_min_kin
                rk_out(k+i,j) = tw%solid_chemistry%rk_old(i)
            end do
            !> Gas kinetic rates not yet implemented in gas_chemistry_c
        else
            !> No kinetic reactions: the reactive mixing reduces to speciating
            !> the conservatively mixed components (u_new == u_tilde).
            u_new(:,j)=u_tilde(:,j)
            call tw%compute_c_nc_from_u_Newton_ideal(tw%get_c1(),u_new(:,j),conc_nc(:,j),niter,CV_flag)
            if (.not.CV_flag) then
                print *, "Target water index: ", tw_idx
                print *, "No convergence in speciation after reactive mixing iteration"
                print *, "Try reducing the time step"
                error stop
            end if
        end if
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
    if (n_kin_tot > 0) then
        write(2,*) "Kinetic reaction rates (rows: reactions [aqueous (linear+redox) | minerals | gas], & 
            first column: reaction name, remaining columns: targets):"
        write(2,*)
        do i=1,n_kin_tot
            write(2,'(2x,a30)',advance='no') react_names(i)(1:min(30,len_trim(react_names(i))))
            write(2,"(*(ES15.5))") (rk_out(i,j), j=1,this%num_target_waters)
        end do
    end if
    close(2)
    deallocate(conc_prior,mix_ratios,u_tilde,u_react,u_new,conc_nc,rk_out,react_names)
contains
    !> Read one full record of arbitrary length from "unit" into an allocatable
    !> string, accumulating fixed-size chunks via non-advancing reads. On normal
    !> end-of-record, ios is returned as 0; on end-of-file or error it is /= 0.
    subroutine read_long_line(unit,line,ios)
        integer(kind=4), intent(in) :: unit
        character(len=:), allocatable, intent(out) :: line
        integer(kind=4), intent(out) :: ios
        character(len=4096) :: chunk
        integer(kind=4) :: sz
        line = ''
        do
            read(unit,'(A)',advance='no',iostat=ios,size=sz) chunk
            if (sz > 0) line = line // chunk(1:sz)
            if (is_iostat_eor(ios)) then
                ios = 0
                exit
            else if (ios /= 0) then
                exit
            end if
        end do
    end subroutine read_long_line

    !> Count whitespace/comma-separated numeric tokens in "line". Treats spaces,
    !> tabs and commas as separators, which is sufficient for list-directed
    !> numeric rows (the mixing-ratio rows).
    integer(kind=4) function count_tokens(line) result(n)
        character(len=*), intent(in) :: line
        integer(kind=4) :: ii, ll
        logical :: in_tok
        character(len=1) :: c
        n = 0
        in_tok = .false.
        ll = len_trim(line)
        do ii = 1, ll
            c = line(ii:ii)
            if (c == ' ' .or. c == ',' .or. c == achar(9)) then
                in_tok = .false.
            else
                if (.not. in_tok) n = n + 1
                in_tok = .true.
            end if
        end do
    end function count_tokens
end subroutine
