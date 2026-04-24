!> \file chem_out_options_m.f90
!> \brief Chemistry output options module.
!> \details
!> Defines the `chem_out_options_s` structure that controls which chemical
!> variables are written to output files during a simulation. The user
!> specifies time steps, water indices, species, reactions, minerals,
!> and variable names to track via an input file (`*_out_opts.dat`).
!>
!> \author Jordi
!> \date Unknown
!> \ingroup io

!> \brief Chemistry output options module.
!> \details Contains the output control structure and its I/O routines.
module chem_out_options_m
    use chem_system_m, only: chem_system_c
    implicit none
    save
    private
    !> \brief Output options structure for chemistry results.
    !> \details Stores user-specified selections of waters, time steps, species,
    !> reactions, minerals, and variable names to include in the output files.
    type, public :: chem_out_options_s
        integer(kind=4) :: num_waters=0                           !< [#] Number of waters to output
        integer(kind=4), allocatable :: ind_waters(:)             !< [-] Indices of waters to output
        integer(kind=4) :: num_time_steps=0                       !< [#] Number of output time steps
        integer(kind=4), allocatable :: time_steps(:)             !< [-] Time step indices for output
        integer(kind=4) :: num_aq_species=0                       !< [#] Number of aqueous species to output
        integer(kind=4), allocatable :: ind_aq_species(:)         !< [-] Indices of aqueous species in phase
        integer(kind=4) :: num_reacts=0                           !< [#] Total number of reactions to output
        integer(kind=4) :: num_reacts_eq=0                        !< [#] Number of equilibrium reactions to output
        integer(kind=4) :: num_reacts_kin=0                       !< [#] Number of kinetic reactions to output
        integer(kind=4), allocatable :: ind_reacts(:)             !< [-] Indices of reactions in chem system (eq then kin)
        integer(kind=4) :: num_mins=0                             !< [#] Number of minerals to output
        integer(kind=4), allocatable :: ind_mins(:)               !< [-] Indices of minerals in chem system
        integer(kind=4) :: num_vars=0                             !< [#] Number of output variable types (-1 = all)
        character(len=256), allocatable :: vars(:)                !< [-] Variable names ('conc', 'react_rate', 'vol_fract')
    contains
        procedure :: set_ind_waters             !< Set water output indices
        procedure :: set_time_steps             !< Set output time step indices
        procedure :: read_chem_out_options      !< Read output options from file
    end type
    
    contains
    
        subroutine set_ind_waters(this,ind_waters)
            implicit none
            class(chem_out_options_s) :: this
            integer(kind=4), intent(in) :: ind_waters(:)
            this%ind_waters=ind_waters
        end subroutine
        
        subroutine set_time_steps(this,time_steps)
            implicit none
            class(chem_out_options_s) :: this
            integer(kind=4), intent(in) :: time_steps(:)
            this%time_steps=time_steps
        end subroutine

        subroutine read_chem_out_options(this,dir,root,chem_syst,num_wat)
            implicit none
            class(chem_out_options_s) :: this
            character(len=*), intent(in) :: dir
            character(len=*), intent(in) :: root
            class(chem_system_c), intent(in) :: chem_syst
            integer(kind=4), intent(in) :: num_wat
            
            integer(kind=4) :: num_time,nits,i,num_mins,num_vars,num_reacts,notw,&
                num_aq_species,flag,j,k
            integer(kind=4) :: unit_local
            integer(kind=4), allocatable :: time_steps(:),ind_wat(:),ind_reacts(:),flags(:)
            character(len=256) :: label,name
            character(len=:), allocatable :: name_trimmed
            
            !num_tar_wat=size(waters)
            num_mins=0
            unit_local = 52
            open(unit=unit_local,file=dir//root//'_out_opts.dat',status='old',action='read')
            do
                read(unit_local,*) label
                if (label=='end') then
                    rewind(unit_local)
                    exit
                else if (label=='TIME STEPS') then
                    read(unit_local,*) nits
                    if (nits<0) then
                        error stop "Invalid number of time steps"
                    else
                        this%num_time_steps=nits+2
                        allocate(time_steps(this%num_time_steps-2))
                        allocate(this%time_steps(this%num_time_steps))
                        if (nits>0) then
                            read(unit_local,*) time_steps
                        end if
                        if (minval(time_steps)<1) then
                            error stop "Invalid time step index"
                        end if
                        this%time_steps(1)=0
                        this%time_steps(2:this%num_time_steps-1)=time_steps
                    end if
                else if (label=='WATERS') then
                    read(unit_local,*) notw
                    if (notw<0) then
                        error stop "Invalid number of waters"
                    else if (notw>0) then
                        this%num_waters=notw
                        allocate(this%ind_waters(this%num_waters))
                        allocate(ind_wat(this%num_waters))
                        read(unit_local,*) ind_wat
                        if (minval(ind_wat)<1 .or. maxval(ind_wat)>num_wat) then
                            error stop "Invalid water index"
                        else
                            this%ind_waters=ind_wat
                        end if
                    else
                        this%num_waters=0
                        !allocate(this%ind_waters(this%num_waters))
                        !do i=1,this%num_waters
                        !    this%ind_waters(i)=i
                        !end do
                    end if
                    !print *, this%ind_waters
                else if (label=='VARIABLES') then
                    read(unit_local,*) num_vars
                    if (num_vars<0) then
                        error stop "Invalid number of variables"
                    else if (num_vars>0) then
                        this%num_vars=num_vars
                        allocate(this%vars(this%num_vars))
                        do i=1,this%num_vars
                            read(unit_local,*) name
                            name_trimmed=trim(name)
                            if (name_trimmed=='conc' .or. name_trimmed=='react_rate' .or. &
                                name_trimmed=='vol_fract') then
                                this%vars(i)=name_trimmed
                            else
                                error stop "Invalid variable name in chem_out_options"
                            end if
                        end do
                    else
                        this%num_vars=0
                    end if
                    !read(unit_local,*) name
                    !if (name=='all') then
                    !    this%num_vars=-1
                    !else if (name=='*') then
                    !    continue
                    !else
                    !    this%num_vars=this%num_vars+1
                    !end if
                else if (label=='AQUEOUS SPECIES') then
                    read(unit_local,*) num_aq_species
                    if (num_aq_species<0 .or. num_aq_species>chem_syst%aq_phase%num_species) then
                        error stop "Invalid number of aqueous species"
                    else if (num_aq_species>0) then
                        this%num_aq_species=num_aq_species
                        allocate(this%ind_aq_species(this%num_aq_species))
                        do i=1,this%num_aq_species
                            read(unit_local,*) name
                            if (trim(name)=='all' .and. this%num_aq_species==num_aq_species) then
                                !this%num_aq_species=num_aq_species
                                forall (i=1:this%num_aq_species) this%ind_aq_species(i)=i
                                exit
                            else if (trim(name)=='all') then
                                error stop "Number of aqueous species does not match 'all' keyword"
                            else    
                                call chem_syst%aq_phase%get_aq_species_index_by_name(trim(name),this%ind_aq_species(i))
                            end if
                        end do
                    else
                        this%num_aq_species=0
                    end if
                else if (label=='MINERALS') then
                    read(unit_local,*) num_mins
                    if (num_mins<0 .or. num_mins>chem_syst%num_minerals) then
                        error stop "Invalid number of minerals"
                    else if (num_mins>0) then
                        this%num_mins=num_mins
                        allocate(this%ind_mins(this%num_mins))
                        do i=1,this%num_mins
                            read(unit_local,*) name
                            if (trim(name)=='all' .and. this%num_mins==num_mins) then
                                !this%num_mins=num_mins
                                forall (i=1:this%num_mins) this%ind_mins(i)=i
                                exit
                            else if (trim(name)=='all') then
                                error stop "Number of minerals does not match 'all' keyword"
                            else
                                call chem_syst%get_mineral_index_by_name(trim(name),this%ind_mins(i))
                            end if
                        end do
                    else
                        this%num_mins=0
                    end if
                else if (label=='REACTIONS') then
                    read(unit_local,*) num_reacts
                    if (num_reacts<0 .or. num_reacts>chem_syst%num_reacts) then
                        error stop "Invalid number of reactions"
                    else if (num_reacts>0) then
                        this%num_reacts=num_reacts
                        allocate(this%ind_reacts(this%num_reacts))
                        allocate(flags(this%num_reacts))
                        allocate(ind_reacts(this%num_reacts))
                        do i=1,this%num_reacts
                            read(unit_local,*) name
                            if (trim(name)=='all' .and. this%num_reacts==num_reacts) then
                                !this%num_reacts=num_reacts
                                this%num_reacts_eq=chem_syst%speciation_alg%num_eq_reactions
                                this%num_reacts_kin=chem_syst%num_kin_reacts
                                forall (i=1:this%num_reacts_eq) this%ind_reacts(i)=i
                                forall (i=1:this%num_reacts_kin) this%ind_reacts(this%num_reacts_eq+i)=i
                                exit
                            else if (trim(name)=='all') then
                                error stop "Number of minerals does not match 'all' keyword"
                            else
                                call chem_syst%is_reaction_in_chem_syst(trim(name),flags(i),ind_reacts(i))
                                if (flags(i)==1) then !> equilibrium reaction found
                                    this%num_reacts_eq=this%num_reacts_eq+1
                                else if (flags(i)==2) then !> kinetic reaction found
                                    this%num_reacts_kin=this%num_reacts_kin+1
                                else
                                    error stop "Reaction name not found in chemical system"
                                end if
                            end if
                        end do
                        j=1 !> counter for equilibrium reactions
                        k=1 !> counter for kinetic reactions
                        do i=1,this%num_reacts
                            if (flags(i)==1) then
                                this%ind_reacts(j)=ind_reacts(i)
                                j=j+1
                            else if (flags(i)==2) then
                                this%ind_reacts(this%num_reacts_eq+k)=ind_reacts(i)
                                k=k+1
                            end if
                            !this%ind_reacts(i)=ind_reacts(i)
                        end do
                    else
                        this%num_reacts=0
                    end if
                else
                    continue
                end if
            end do
            close(unit_local)
        end subroutine
end module