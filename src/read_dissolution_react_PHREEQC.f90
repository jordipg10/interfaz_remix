subroutine read_dissolution_react_PHREEQC(this,string,phase)
    use eq_reaction_m, only: eq_reaction_c
    use phase_m, only: phase_c
    use species_m, only: species_c
    use strings_m
    implicit none
    class(eq_reaction_c) :: this
    character(len=*), intent(in) :: string !> dissolution reaction string
    class(phase_c), intent(inout), optional :: phase !> defined phase
    
    integer(kind=4) :: i,j,k,n,num_reactants,n_ions,eq_ind,ind_init,plus_ind,minus_ind,sign_ind,phase_ind
    real(kind=8) :: logK
    real(kind=8), allocatable :: electr_valences(:)
    character(len=256) :: str
    character(len=256), allocatable :: str_arr(:),species_names(:)
    logical :: gas_flag
    type(species_c), allocatable :: local_species(:)
    
    
    n=len_trim(string)

    ind_init=1
    do i=1,n
        if (string(i:i)=='') then
            ind_init=ind_init+1
        else
            exit
        end if
    end do
    this%num_species=1
    eq_ind=index(string,'=')
    
    do i=ind_init,n-2
        if (string(i:i+2)==' + ' .or. string(i:i+2)==' - ' .or. string(i:i+2)==' = ') then
            this%num_species=this%num_species+1
        else
            continue
        end if
    end do
    call this%allocate_reaction()
    allocate(local_species(this%num_species))
    
    plus_ind=index(string(ind_init:n),'+')
    minus_ind=index(string(ind_init:n),'-')
    if (plus_ind>0 .and. minus_ind>0) then
        sign_ind=min(plus_ind,minus_ind)
    else if (plus_ind>0) then
        sign_ind=plus_ind
    else if (minus_ind>0) then
        sign_ind=minus_ind
    else
        sign_ind=eq_ind
    end if
    phase_ind=min(sign_ind,eq_ind)

    if (present(phase)) then
        call phase%set_phase_name(string(ind_init:phase_ind-2))
    end if
    
    j=1
    k=ind_init
    i=ind_init
    do while (j<this%num_species)
        if (string(i:i+2)==' + ' .or. string(i:i+2)==' - ' .or. string(i:i+2)==' = ') then
            call local_species(j)%set_name(string(k:i-1))
            if (string(i-1:i-1)=='+') then
                call local_species(j)%set_valence(1)
            else if (string(i-1:i-1)=='-') then
                call local_species(j)%set_valence(-1)
            else if (string(i-2:i-2)=='+' .or. string(i-2:i-2)=='-') then
                read(string(i-1:i-1),*) local_species(j)%valence
            end if
            k=i+3
            j=j+1
        end if
        i=i+1
    end do
    call local_species(j)%set_name(string(k:n))
    if (string(n-1:n-1)=='+') then
        call local_species(j)%set_valence(1)
    else if (string(n-1:n-1)=='-') then
        call local_species(j)%set_valence(-1)
    else if (string(n-2:n-2)=='+' .or. string(n-2:n-2)=='-') then
        read(string(n-1:n-1),*) local_species(j)%valence
    end if

end subroutine