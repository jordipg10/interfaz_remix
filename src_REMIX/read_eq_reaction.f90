subroutine read_eq_reaction(this,species,filename)
    use eq_reaction_m, only: eq_reaction_c, species_c
    use strings_m
    implicit none
    class(eq_reaction_c) :: this
    class(species_c), intent(in) :: species !> defines equilibrium reaction
    character(len=*), intent(in) :: filename
    
    integer(kind=4) :: i,j,n_react, n_ions
    real(kind=8) :: logK
    real(kind=8), allocatable :: electr_valences(:)
    character(len=256) :: str
    character(len=256), allocatable :: str_arr(:),species_names(:)
    logical :: flag
    

    open(unit=25,file=filename,status='old',action='read')
!> Archivo ad hoc
    if (filename=='reacciones_disociacion.dat') then
        do
            read(25,*) str, n_react
            if (str==species%name) then
                call this%allocate_reaction(n_react+1)
                backspace(25)
                read(25,*) str, n_react, (this%stoichiometry(j), this%species(j)%name, j=1,n_react), logK
                exit
            else if (str=='end') then
                exit
            else
                continue
            end if
        end do
        do j=1,n_react
            this%species(j)%cst_act_flag=.false.
            this%stoichiometry(j)=-this%stoichiometry(j)
        end do
        call this%species(n_react+1)%assign_species(species)
        this%stoichiometry(n_react+1)=1d0
        this%eq_cst=10**(logK)
    !else if (filename=='reacciones_adsorcion.dat') then
    !>    n_ions=size(species)
    !>    allocate(str_arr(n_ions))
    !>    do
    !>        read(25,*) str
    !>        if (str==species(1)%name) then
    !>            backspace(25)
    !>            read(25,*) (str_arr(j), j=1,n_ions)
    !>            call compare_str_arrays(str_arr,species_names,flag)
    !>            if (flag==.true.) then
    !>                call this%allocate_reaction(2*n_ions)
    !>                allocate(electr_valences(n_ions))
    !>                backspace(25)
    !>                read(25,*) (str_arr(j), j=1,n_ions), (electr_valences(j), j=1,n_ions), this%eq_cst
    !>                exit
    !>            else
    !>                continue
    !>            end if
    !>        else if (str=='end') then
    !>            exit
    !>        else
    !>            continue
    !>        end if
    !>    end do
    !>    do j=1,n_ions
    !>        this%stoichiometry(j)=-1d0/electr_valences(j)
    !>        this%species(j)%name=species(j)%name
    !>        this%stoichiometry(n_ions+1)=1d0/electr_valences(j)
    !>        this%species(n_ions+1)%name='X-'//species(j)%name
    !>    end do
    !>    !this%stoichiometry(n_react+1)%name=str
    !>    !this%stoichiometry(n_react+1)%stoich_coeff=1d0
    !>    !this%eq_cst=10**(logK)
    end if
    close(25)
end subroutine