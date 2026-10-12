!> @file read_PHREEQC_DB.f90
!> @brief Reads a PHREEQC thermodynamic database (e.g. 'phreeqc.dat') into an
!>        already-defined chemical system.
!>
!> @details
!> This subroutine assumes that the SPECIES of the chemical system are ALREADY known
!> before the database is read, i.e. the following must be populated (with names) on
!> entry:
!>   - this%species           : global species array (aqueous, minerals, gases, surface)
!>   - this%aq_phase%aq_species
!>   - this%minerals
!>   - this%gas_phase%gases
!>   - this%eq_reacts         : equilibrium reactions, each with %name set to the name
!>                              of the species / mineral / gas / exchange complex it defines
!>
!> The reader then scans the PHREEQC database block by block and, for every entry that
!> belongs to the chemical system (matched BY NAME), fills:
!>   - eq_cst_ref     from '-log_k'      (stored as K = 10**log_k)
!>   - delta_h_ref    from '-delta_h'    (converted to J/mol; PHREEQC default kJ/mol, 'kcal'/'cal'/'J' honoured)
!>   - coeffs_logK_T  from '-analytic' / '-analytical_expression' (6 coeffs, zero padded)
!>   - reaction stoichiometry + species indices (when all species are resolvable)
!>   - aqueous species activity parameters from '-gamma' and Redlich volume from '-Vm'
!>   - mineral molar volume from '-Vm' (cm3/mol -> L/mol)
!>   - gas critical properties from '-T_c', '-P_c', '-Omega'
!>   - primary aqueous species valence from SOLUTION_MASTER_SPECIES
!>
!> Blocks parsed: SOLUTION_MASTER_SPECIES, SOLUTION_SPECIES, PHASES, EXCHANGE_SPECIES,
!> SURFACE_SPECIES. Other blocks (EXCHANGE_MASTER_SPECIES, SURFACE_MASTER_SPECIES, RATES,
!> END, ...) are skipped. '-dw' (diffusion) is parsed over but not stored (no field).
!>
!> Sign convention: the DEFINED species (secondary aqueous species / mineral / gas /
!> exchange complex) is given a NEGATIVE stoichiometric coefficient and the species it is
!> made of POSITIVE coefficients, matching the rest of REMIX (see read_master25).
!>
!> @author Jordi Petchamé-Guerrero
!> @date 2026
subroutine read_PHREEQC_DB(this,filename)
    use chem_system_m, only: chem_system_c
    implicit none
    class(chem_system_c) :: this          !< Chemical system (species already known)
    character(len=*), intent(in) :: filename  !< Full path to the PHREEQC database file

    !> Maximum number of species allowed in a single reaction
    integer(kind=4), parameter :: MAXSP=40
    real(kind=8), parameter :: CAL_TO_J=4.184d0       !< cal -> J conversion factor

    integer(kind=4) :: unit,ios,block_id,lbl,i
    character(len=512) :: raw
    character(len=:), allocatable :: line,sline
    character(len=64) :: tok1

    !> --- current-entry state (host associated with the internal procedures) ---
    integer(kind=4) :: nsp,rhs_start,cur_block
    character(len=64) :: rnames(MAXSP)
    real(kind=8) :: rstoich(MAXSP)
    character(len=64) :: cur_name
    logical :: have_entry,cur_has_react,cur_is_gas
    !> parsed thermodynamic / property data for the current entry
    real(kind=8) :: e_logk,e_dh,e_coeffs(6),e_gamma_a,e_gamma_b,e_vm(10),e_Tc,e_Pc,e_Om
    integer(kind=4) :: e_nvm
    logical :: e_has_logk,e_has_dh,e_has_analytic,e_has_gamma,e_has_vm
    logical :: e_has_Tc,e_has_Pc,e_has_Om

    !> Pre-condition: the species must be known before reading the database
    if (.not. allocated(this%species)) then
        error stop "read_PHREEQC_DB: chemical system species must be defined before reading the database"
    end if

    open(newunit=unit,file=filename,form='formatted',status='old',action='read',iostat=ios)
    if (ios/=0) then
        error stop "read_PHREEQC_DB: could not open database file '"//trim(filename)//"'"
    end if

    block_id=0        !< 0=none, 1=SMS, 2=SS, 3=PHASES, 4=EXCH, 5=SURF, 9=skip
    have_entry=.false.
    call reset_entry()

    do
        read(unit,'(A)',iostat=ios) raw
        if (ios/=0) exit                    !< end of file
        line=decomment(raw)                 !< drop trailing '#' comment
        sline=lstrip(line)                  !< remove leading spaces/tabs
        if (len_trim(sline)==0) cycle       !< blank / comment-only line
        tok1=first_token(sline)

        !> Block label ?
        lbl=block_label_id(tok1)
        if (lbl/=0) then
            call flush_entry()
            block_id=lbl
            cycle
        end if

        select case (block_id)
        case (1)   !> SOLUTION_MASTER_SPECIES (line-based)
            call process_master_species(sline)
        case (2,4,5)   !> SOLUTION_SPECIES / EXCHANGE_SPECIES / SURFACE_SPECIES
            if (index(sline,'=')>0) then
                call flush_entry()
                call start_entry(block_id)
                call parse_reaction(sline)
                cur_has_react=.true.
                if (nsp>=1) cur_name=rnames(rhs_start)   !< defined = first product
            else
                call handle_property_line(sline)
            end if
        case (3)   !> PHASES (name header, then reaction, then properties)
            if (index(sline,'=')>0) then
                call parse_reaction(sline)               !< reaction of current phase
                cur_has_react=.true.
            else if (is_prop_line(sline)) then
                call handle_property_line(sline)
            else
                call flush_entry()
                call start_entry(3)
                cur_name=tok1                            !< phase name (defines the reaction)
                cur_is_gas=(index(cur_name,'(g)')>0)
            end if
        case default
            continue                                     !< inside a skipped block
        end select
    end do
    call flush_entry()                                   !< flush the last entry
    close(unit)

contains

    !> Reset the current-entry accumulators
    subroutine reset_entry()
        nsp=0; rhs_start=1; cur_name=''; cur_block=0
        cur_has_react=.false.; cur_is_gas=.false.
        e_logk=0d0; e_dh=0d0; e_coeffs=0d0
        e_gamma_a=0d0; e_gamma_b=0d0; e_vm=0d0; e_nvm=0
        e_Tc=0d0; e_Pc=0d0; e_Om=0d0
        e_has_logk=.false.; e_has_dh=.false.; e_has_analytic=.false.
        e_has_gamma=.false.; e_has_vm=.false.
        e_has_Tc=.false.; e_has_Pc=.false.; e_has_Om=.false.
    end subroutine

    !> Begin a new entry in the given block
    subroutine start_entry(blk)
        integer(kind=4), intent(in) :: blk
        call reset_entry()
        cur_block=blk
        have_entry=.true.
    end subroutine

    !> Apply the parsed current entry to the chemical system
    subroutine flush_entry()
        if (.not. have_entry) return
        select case (cur_block)
        case (2)   !> aqueous species / complex
            if (cur_has_react) then
                if (is_primary_reaction()) then
                    call fill_aq_species_props(cur_name)         !< "X = X": primary species
                else
                    call fill_aq_species_props(cur_name)         !< the complex
                    call fill_reaction(1)                        !< aqueous complexation
                end if
            end if
        case (3)   !> PHASES
            if (cur_is_gas) then
                call fill_gas_props(cur_name)
                call fill_reaction(6)                            !< gas dissolution
            else
                call fill_mineral_props(cur_name)
                call fill_reaction(2)                            !< mineral dissolution
            end if
        case (4,5) !> exchange / surface complexation
            if (cur_has_react .and. .not. is_primary_reaction()) then
                call fill_reaction(3)                            !< (surface) complexation
            end if
        end select
        have_entry=.false.
    end subroutine

    !> TRUE if the current reaction is trivial ("X = X"): a primary/master species
    logical function is_primary_reaction()
        is_primary_reaction=(nsp==2 .and. trim(rnames(1))==trim(rnames(2)))
    end function

    !> Fill an existing equilibrium reaction (matched by defined-species name)
    subroutine fill_reaction(rtype)
        integer(kind=4), intent(in) :: rtype
        integer(kind=4) :: ridx,i,p,idxs(MAXSP)
        logical :: allfound
        character(len=64) :: lname

        ridx=find_eq_react(cur_name)
        if (ridx==0) return                  !< reaction not in this chemical system

        !> resolve every reaction species to a global species index
        allfound=.true.
        do i=1,nsp
            lname=rnames(i)
            if (cur_block==3 .and. i==1) lname=cur_name   !< PHASES: 1st term is the phase
            idxs(i)=find_global_species(lname)
            if (idxs(i)==0) allfound=.false.
        end do

        if (allfound) then
            !> position of the defined species (negative coefficient)
            if (cur_block==3) then
                p=1
            else
                p=rhs_start
            end if
            if (rstoich(p)>0d0) rstoich(1:nsp)=-rstoich(1:nsp)
            if (allocated(this%eq_reacts(ridx)%species_ind)) then
                call this%eq_reacts(ridx)%deallocate_reaction()
            end if
            call this%eq_reacts(ridx)%allocate_reaction(nsp)
            do i=1,nsp
                call this%eq_reacts(ridx)%set_single_species(i,idxs(i))
            end do
            call this%eq_reacts(ridx)%set_stoichiometry(rstoich(1:nsp))
            call this%eq_reacts(ridx)%set_react_type(rtype)
            call this%eq_reacts(ridx)%set_react_name(cur_name)
        end if

        !> thermodynamic data (independent of stoichiometry resolution)
        if (e_has_logk) call this%eq_reacts(ridx)%set_eq_cst_ref(10d0**e_logk)
        if (e_has_dh)   call this%eq_reacts(ridx)%set_delta_h_ref(e_dh)
        if (e_has_analytic) then
            if (.not. allocated(this%eq_reacts(ridx)%coeffs_logK_T)) then
                allocate(this%eq_reacts(ridx)%coeffs_logK_T(6))
            end if
            this%eq_reacts(ridx)%coeffs_logK_T=e_coeffs
        end if
    end subroutine

    !> Fill activity ('-gamma') and Redlich volume ('-Vm') parameters of an aqueous species
    subroutine fill_aq_species_props(nm)
        character(len=*), intent(in) :: nm
        integer(kind=4) :: aidx,gidx
        aidx=find_aq_species(nm)
        if (aidx>0) then
            if (e_has_gamma) then
                this%aq_phase%aq_species(aidx)%params_act_coeff%ion_size_param=e_gamma_a
                this%aq_phase%aq_species(aidx)%params_act_coeff%a_TJ=e_gamma_a
                this%aq_phase%aq_species(aidx)%params_act_coeff%b_TJ=e_gamma_b
            end if
            if (e_has_vm) then
                call this%aq_phase%aq_species(aidx)%params_spec_vol_Redlich%set_params(e_vm)
            end if
        end if
        gidx=find_global_species(nm)
        if (gidx>0) then
            if (e_has_gamma) then
                this%species(gidx)%params_act_coeff%ion_size_param=e_gamma_a
                this%species(gidx)%params_act_coeff%a_TJ=e_gamma_a
                this%species(gidx)%params_act_coeff%b_TJ=e_gamma_b
            end if
            if (e_has_vm) then
                call this%species(gidx)%params_spec_vol_Redlich%set_params(e_vm)
            end if
        end if
    end subroutine

    !> Fill molar volume of a mineral from '-Vm' (cm3/mol -> dm3/mol)
    subroutine fill_mineral_props(nm)
        character(len=*), intent(in) :: nm
        integer(kind=4) :: midx,gidx
        if (.not. e_has_vm) return
        midx=find_mineral(nm)
        if (midx>0) call this%minerals(midx)%mineral%set_mol_vol(e_vm(1)*1d-3)
        gidx=find_global_species(nm)
        if (gidx>0) call this%species(gidx)%set_mol_vol(e_vm(1)*1d-3)
    end subroutine

    !> Fill critical properties of a gas from '-T_c', '-P_c', '-Omega'
    subroutine fill_gas_props(nm)
        character(len=*), intent(in) :: nm
        integer(kind=4) :: gidx
        gidx=find_gas(nm)
        if (gidx==0) return
        if (e_has_Tc) call this%gas_phase%gases(gidx)%set_crit_temp(e_Tc)
        if (e_has_Pc) call this%gas_phase%gases(gidx)%set_crit_press(e_Pc)
        if (e_has_Om) call this%gas_phase%gases(gidx)%set_acentric_fact(e_Om)
    end subroutine

    !> Set the valence of a primary aqueous species from SOLUTION_MASTER_SPECIES
    subroutine process_master_species(s)
        character(len=*), intent(in) :: s
        character(len=64) :: sp_tok
        integer(kind=4) :: aidx,gidx,v
        sp_tok=first_token(after_key(s))     !< second column = master species
        if (len_trim(sp_tok)==0) return
        v=parse_charge(sp_tok)
        aidx=find_aq_species(sp_tok)
        if (aidx>0) call this%aq_phase%aq_species(aidx)%set_valence(v)
        gidx=find_global_species(sp_tok)
        if (gidx>0) call this%species(gidx)%set_valence(v)
    end subroutine

    !> Split a property line on ';' and process each keyword segment
    subroutine handle_property_line(s)
        character(len=*), intent(in) :: s
        integer(kind=4) :: i,n,start
        n=len_trim(s)
        start=1
        do i=1,n
            if (s(i:i)==';') then
                call process_property(lstrip(s(start:i-1)))
                start=i+1
            end if
        end do
        call process_property(lstrip(s(start:n)))
    end subroutine

    !> Parse one property keyword segment into the current-entry accumulators
    subroutine process_property(seg)
        character(len=*), intent(in) :: seg
        character(len=64) :: key
        character(len=:), allocatable :: k,rest
        real(kind=8) :: tmp(10)
        integer(kind=4) :: n,ios2
        if (len_trim(seg)==0) return
        key=first_token(seg)
        k=to_lower(trim(key))
        if (len(k)>=1) then
            if (k(1:1)=='-') k=k(2:len(k))        !< strip leading '-'
        end if
        rest=after_key(seg)
        select case (k)
        case ('log_k','logk')
            read(rest,*,iostat=ios2) e_logk
            if (ios2==0) e_has_logk=.true.
        case ('delta_h')
            call get_reals(rest,tmp,1,n)
            if (n>=1) then
                !> stored in J/mol; the unit token (if any) follows the value
                block
                    character(len=:), allocatable :: u
                    u=to_lower(after_key(rest))
                    if (index(u,'kcal')>0) then
                        e_dh=tmp(1)*CAL_TO_J*1d3
                    else if (index(u,'cal')>0) then
                        e_dh=tmp(1)*CAL_TO_J
                    else if (index(u,'kj')>0) then
                        e_dh=tmp(1)*1d3
                    else if (index(u,'j')>0) then
                        e_dh=tmp(1)
                    else
                        e_dh=tmp(1)*1d3                !< kJ/mol (PHREEQC default)
                    end if
                end block
                e_has_dh=.true.
            end if
        case ('analytic','analytical','analytical_expression')
            call get_reals(rest,tmp,6,n)
            e_coeffs=0d0
            if (n>=1) e_coeffs(1:n)=tmp(1:n)
            e_has_analytic=.true.
        case ('gamma')
            call get_reals(rest,tmp,2,n)
            if (n>=1) e_gamma_a=tmp(1)
            if (n>=2) e_gamma_b=tmp(2)
            e_has_gamma=.true.
        case ('vm')
            call get_reals(rest,e_vm,10,e_nvm)
            if (e_nvm>=1) e_has_vm=.true.
        case ('t_c')
            read(rest,*,iostat=ios2) e_Tc
            if (ios2==0) e_has_Tc=.true.
        case ('p_c')
            read(rest,*,iostat=ios2) e_Pc
            if (ios2==0) e_has_Pc=.true.
        case ('omega')
            read(rest,*,iostat=ios2) e_Om
            if (ios2==0) e_has_Om=.true.
        case default
            continue                               !< '-dw' and unknown keywords ignored
        end select
    end subroutine

    !> Parse a reaction string "LHS = RHS" into rnames/rstoich (LHS<0, RHS>0)
    subroutine parse_reaction(s)
        character(len=*), intent(in) :: s
        integer(kind=4) :: eqpos
        eqpos=index(s,'=')
        nsp=0
        call parse_side(s(1:eqpos-1),-1d0)
        rhs_start=nsp+1
        call parse_side(s(eqpos+1:len_trim(s)),1d0)
    end subroutine

    !> Split one side of a reaction on ' + ' and append its terms
    subroutine parse_side(side,sgn)
        character(len=*), intent(in) :: side
        real(kind=8), intent(in) :: sgn
        integer(kind=4) :: n,i,start,sep
        real(kind=8) :: coeff
        character(len=64) :: nm
        n=len_trim(side)
        if (n==0) return
        start=1
        do
            sep=0
            do i=start,n-2
                if (side(i:i+2)==' + ') then
                    sep=i
                    exit
                end if
            end do
            if (sep==0) then
                call parse_term(side(start:n),coeff,nm)
                if (len_trim(nm)>0 .and. nsp<MAXSP) then
                    nsp=nsp+1; rnames(nsp)=nm; rstoich(nsp)=sgn*coeff
                end if
                exit
            else
                call parse_term(side(start:sep-1),coeff,nm)
                if (len_trim(nm)>0 .and. nsp<MAXSP) then
                    nsp=nsp+1; rnames(nsp)=nm; rstoich(nsp)=sgn*coeff
                end if
                start=sep+3
            end if
        end do
    end subroutine

    !> Parse a single reaction term "[coeff ]name" (coeff may be attached, e.g. "2H2O")
    subroutine parse_term(term,coeff,nm)
        character(len=*), intent(in) :: term
        real(kind=8), intent(out) :: coeff
        character(len=*), intent(out) :: nm
        character(len=:), allocatable :: t
        integer(kind=4) :: L,j,ios2
        t=lstrip(term)
        L=len_trim(t)
        coeff=1d0; nm=''
        if (L==0) return
        j=1
        do while (j<=L .and. ((t(j:j)>='0' .and. t(j:j)<='9') .or. t(j:j)=='.'))
            j=j+1
        end do
        if (j>1) then
            read(t(1:j-1),*,iostat=ios2) coeff
            if (ios2/=0) coeff=1d0
        end if
        do while (j<=L .and. t(j:j)==' ')
            j=j+1
        end do
        if (j<=L) nm=t(j:L)
    end subroutine

    !> Read up to maxn real numbers from a whitespace-separated string
    subroutine get_reals(s,vals,maxn,nout)
        character(len=*), intent(in) :: s
        real(kind=8), intent(out) :: vals(:)
        integer(kind=4), intent(in) :: maxn
        integer(kind=4), intent(out) :: nout
        integer(kind=4) :: i,j,L,ios2
        real(kind=8) :: rval
        nout=0
        L=len_trim(s)
        i=1
        do while (i<=L .and. nout<maxn)
            do while (i<=L .and. (s(i:i)==' ' .or. s(i:i)==char(9)))
                i=i+1
            end do
            if (i>L) exit
            j=i
            do while (j<=L .and. .not.(s(j:j)==' ' .or. s(j:j)==char(9)))
                j=j+1
            end do
            read(s(i:j-1),*,iostat=ios2) rval
            if (ios2/=0) exit
            nout=nout+1
            vals(nout)=rval
            i=j
        end do
    end subroutine

    !> Extract the integer charge encoded at the end of a species name (e.g. "CO3-2" -> -2)
    integer(kind=4) function parse_charge(nm) result(v)
        character(len=*), intent(in) :: nm
        integer(kind=4) :: L,j,mag,ios2
        v=0
        L=len_trim(nm)
        if (L==0) return
        if (nm(L:L)>='0' .and. nm(L:L)<='9') then
            j=L
            do while (j>=1 .and. nm(j:j)>='0' .and. nm(j:j)<='9')
                j=j-1
            end do
            if (j>=1) then
                read(nm(j+1:L),*,iostat=ios2) mag
                if (ios2/=0) mag=0
                if (nm(j:j)=='+') then
                    v=mag
                else if (nm(j:j)=='-') then
                    v=-mag
                end if
            end if
        else if (nm(L:L)=='+') then
            v=1
        else if (nm(L:L)=='-') then
            v=-1
        end if
    end function

    !> --- name-based lookups (return 0 if not found) ---

    integer(kind=4) function find_global_species(nm) result(idx)
        character(len=*), intent(in) :: nm
        integer(kind=4) :: i
        idx=0
        do i=1,size(this%species)
            if (trim(this%species(i)%name)==trim(nm)) then
                idx=i; return
            end if
        end do
    end function

    integer(kind=4) function find_aq_species(nm) result(idx)
        character(len=*), intent(in) :: nm
        integer(kind=4) :: i
        idx=0
        if (.not. allocated(this%aq_phase%aq_species)) return
        do i=1,this%aq_phase%num_species
            if (trim(this%aq_phase%aq_species(i)%name)==trim(nm)) then
                idx=i; return
            end if
        end do
    end function

    integer(kind=4) function find_mineral(nm) result(idx)
        character(len=*), intent(in) :: nm
        integer(kind=4) :: i
        idx=0
        do i=1,this%num_minerals
            if (trim(this%minerals(i)%mineral%name)==trim(nm) .or. &
                trim(this%minerals(i)%name)==trim(nm)) then
                idx=i; return
            end if
        end do
    end function

    integer(kind=4) function find_gas(nm) result(idx)
        character(len=*), intent(in) :: nm
        integer(kind=4) :: i
        idx=0
        if (.not. allocated(this%gas_phase%gases)) return
        do i=1,size(this%gas_phase%gases)
            if (trim(this%gas_phase%gases(i)%name)==trim(nm)) then
                idx=i; return
            end if
        end do
    end function

    integer(kind=4) function find_eq_react(nm) result(idx)
        character(len=*), intent(in) :: nm
        integer(kind=4) :: i
        idx=0
        if (.not. allocated(this%eq_reacts)) return
        do i=1,size(this%eq_reacts)
            if (trim(this%eq_reacts(i)%name)==trim(nm)) then
                idx=i; return
            end if
        end do
    end function

    !> --- string helpers ---

    !> Return the PHREEQC block id of a token (0 = not a block label, 9 = skipped block)
    integer(kind=4) function block_label_id(tok) result(id)
        character(len=*), intent(in) :: tok
        character(len=:), allocatable :: t
        t=trim(tok)
        select case (t)
        case ('SOLUTION_MASTER_SPECIES'); id=1
        case ('SOLUTION_SPECIES');        id=2
        case ('PHASES');                  id=3
        case ('EXCHANGE_SPECIES');        id=4
        case ('SURFACE_SPECIES');         id=5
        case ('EXCHANGE_MASTER_SPECIES','SURFACE_MASTER_SPECIES','RATES','END', &
              'KNOBS','PITZER','SIT','LLNL_AQUEOUS_MODEL_PARAMETERS','KINETICS')
            id=9
        case default
            id=0
        end select
    end function

    !> TRUE if a line is a property line ('-keyword ...' or a bare recognised keyword)
    logical function is_prop_line(s)
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: k
        is_prop_line=.false.
        if (len_trim(s)==0) return
        if (s(1:1)=='-') then
            is_prop_line=.true.; return
        end if
        k=to_lower(trim(first_token(s)))
        select case (k)
        case ('log_k','logk','delta_h','analytic','analytical','analytical_expression', &
              'gamma','dw','vm','t_c','p_c','omega')
            is_prop_line=.true.
        end select
    end function

    !> Remove any trailing '#'-comment from a line
    function decomment(s) result(r)
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: r
        integer(kind=4) :: h
        h=index(s,'#')
        if (h>0) then
            r=s(1:h-1)
        else
            r=s
        end if
    end function

    !> Remove leading spaces and tabs
    function lstrip(s) result(r)
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: r
        integer(kind=4) :: i,n
        n=len(s)
        i=1
        do while (i<=n .and. (s(i:i)==' ' .or. s(i:i)==char(9)))
            i=i+1
        end do
        if (i>len_trim(s)) then
            r=''
        else
            r=s(i:len_trim(s))
        end if
    end function

    !> First whitespace-delimited token of a string
    function first_token(s) result(tok)
        character(len=*), intent(in) :: s
        character(len=64) :: tok
        character(len=:), allocatable :: t
        integer(kind=4) :: i,L
        t=lstrip(s)
        L=len_trim(t)
        tok=''
        i=1
        do while (i<=L .and. .not.(t(i:i)==' ' .or. t(i:i)==char(9)))
            i=i+1
        end do
        if (i>1) tok=t(1:i-1)
    end function

    !> Everything after the first token (trimmed on the left)
    function after_key(s) result(r)
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: r,t
        integer(kind=4) :: i,L
        t=lstrip(s)
        L=len_trim(t)
        i=1
        do while (i<=L .and. .not.(t(i:i)==' ' .or. t(i:i)==char(9)))
            i=i+1
        end do
        if (i>L) then
            r=''
        else
            r=lstrip(t(i:L))
        end if
    end function

    !> Lowercase a string
    function to_lower(s) result(r)
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: r
        integer(kind=4) :: i,c
        r=s
        do i=1,len(r)
            c=iachar(r(i:i))
            if (c>=iachar('A') .and. c<=iachar('Z')) r(i:i)=achar(c+32)
        end do
    end function

end subroutine
