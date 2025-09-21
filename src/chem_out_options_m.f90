!> Chemical output type: contains variables that will be considered in output file
module chem_out_options_m
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
    save
    type, public :: chem_out_options_t
        integer(kind=4) :: num_target_waters
        integer(kind=4), allocatable :: ind_target_waters(:)
        integer(kind=4) :: num_time_steps
        integer(kind=4), allocatable :: time_steps(:)
        integer(kind=4) :: num_aq_species
        integer(kind=4), allocatable :: ind_aq_species(:)
        integer(kind=4), allocatable :: num_reacts(:)
        integer(kind=4), allocatable :: ind_reacts(:)
        integer(kind=4), allocatable :: num_mins(:)
        integer(kind=4), allocatable :: ind_mins(:)
        integer(kind=4) :: num_vars !> if -1: all variables will be written (not recommended)
        character(len=256), allocatable :: vars(:) !> variables
    contains
        procedure, private :: set_ind_target_waters
        procedure, private :: set_time_steps
        procedure, public :: read_chem_out_options
    end type
    
    contains
    
        subroutine set_ind_target_waters(this,ind_target_waters)
            implicit none
            class(chem_out_options_t) :: this
            integer(kind=4), intent(in) :: ind_target_waters(:)
            this%ind_target_waters=ind_target_waters
        end subroutine
        
        subroutine set_time_steps(this,time_steps)
            implicit none
            class(chem_out_options_t) :: this
            integer(kind=4), intent(in) :: time_steps(:)
            this%time_steps=time_steps
        end subroutine
        
        subroutine read_chem_out_options(this,root,target_waters)
            implicit none
            class(chem_out_options_t) :: this
            character(len=*), intent(in) :: root
            !integer(kind=4), intent(in) :: unit
            type(aqueous_chemistry_c), intent(in) :: target_waters(:)
            
            integer(kind=4) :: num_time,nits,i,num_mins,num_vars,num_reacts,unit
            character(len=256) :: label,name
            
            this%num_aq_species=0
            num_mins=0
            this%num_vars=0
            open(unit=5,file=root//'_out_options.dat',status='old',action='read')
            do
                read(5,*) label
                if (label=='end') then
                    rewind(5)
                    exit
                else if (label=='TIME STEPS') then
                    read(5,*) nits
                    if (nits<0) then
                        error stop "Invalid number of time steps"
                    else
                        this%num_time_steps=nits+2
                        allocate(this%time_steps(this%num_time_steps))
                        if (nits>0) then
                            read(5,*) this%time_steps(2:this%num_time_steps-1)
                        end if
                        this%time_steps(1)=0
                    end if
                    !write(*,*) this%time_steps
                else if (label=='TARGET WATERS') then
                    read(5,*) this%num_target_waters
                    if (this%num_target_waters<0) then
                        error stop "Invalid number of target waters"
                    else if (this%num_target_waters>0) then
                        allocate(this%ind_target_waters(this%num_target_waters))
                        read(5,*) this%ind_target_waters
                    else
                        this%num_target_waters=size(target_waters)
                        allocate(this%ind_target_waters(this%num_target_waters))
                        do i=1,this%num_target_waters
                            this%ind_target_waters(i)=i
                        end do
                    end if
                    !print *, this%ind_target_waters
                else if (label=='VARIABLES') then
                    read(5,*) name
                    if (name=='all') then
                        this%num_vars=-1
                    else if (name=='*') then
                        continue
                    else
                        this%num_vars=this%num_vars+1
                    end if
                else if (label=='AQUEOUS VARIABLE ACTIVITY SPECIES') then
                    read(5,*) name
                    if (name=='all') then
                        this%num_aq_species=target_waters(1)%aq_phase%num_species-target_waters(1)%aq_phase%wat_flag
                        allocate(this%ind_aq_species(this%num_aq_species))
                        do i=1,this%num_aq_species
                            this%ind_aq_species(i)=target_waters(1)%aq_phase%ind_diss_solids(i)
                        end do
                    else if (name=='*') then
                        continue
                    else
                        this%num_aq_species=this%num_aq_species+1
                    end if
                else if (label=='MINERALS') then
                    allocate(this%num_mins(this%num_target_waters))
                    read(5,*) name
                    if (name=='all') then
                        do i=1,this%num_target_waters
                            this%num_mins(i)=target_waters(this%ind_target_waters(i))%solid_chemistry%reactive_zone%num_minerals
                        end do
                    else if (name=='*') then
                        continue
                    else
                        num_mins=num_mins+1
                    end if
                else if (label=='REACTIONS') then
                    allocate(this%num_reacts(this%num_target_waters))
                    read(5,*) name
                    if (name=='all') then
                        do i=1,this%num_target_waters
                            this%num_reacts(i)=target_waters(&
                            this%ind_target_waters(i))%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions + &
                            target_waters(this%ind_target_waters(i))%solid_chemistry%reactive_zone%chem_syst%num_kin_reacts
                        end do
                    else if (name=='*') then
                        continue
                    else
                        num_reacts=num_reacts+1
                    end if
                else
                    continue
                end if
            end do
            close(5)
        end subroutine
end module