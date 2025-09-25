subroutine read_association_react_PHREEQC(this,string,prim_flag,defined_species)
    use eq_reaction_m, only : eq_reaction_c, aq_species_c
    implicit none
    class(eq_reaction_c) :: this
    character(len=*), intent(in) :: string !> association reaction
    logical, intent(out) :: prim_flag !> TRUE if species is primary
    type(aq_species_c), intent(out), optional :: defined_species !> defined aqueous species
    
    integer(kind=4) :: i,j,k,n,num_reactants,n_ions,eq_ind,ind_init,plus_ind,minus_ind,sign_ind,hash_ind
    real(kind=8) :: logK
    real(kind=8), allocatable :: electr_valences(:)
    character(len=256) :: str
    character(len=256), allocatable :: str_arr(:),species_names(:)
    logical :: flag
    
    
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
    hash_ind=index(string,'#')
    if (string(ind_init:eq_ind-2)==string(eq_ind+2:n) .and. hash_ind==0) then
        prim_flag=.true.
        if (present(defined_species)) then
            call defined_species%set_name(string(ind_init:eq_ind-2))
        end if
    else if (string(ind_init:eq_ind-2)==string(eq_ind+2:hash_ind-2) .and. hash_ind>0) then
        prim_flag=.true.
        if (present(defined_species)) then
            call defined_species%set_name(string(ind_init:eq_ind-2))
        end if
    else
        prim_flag=.false.
        do i=ind_init,n-2
            if (string(i:i+2)==' + ' .or. string(i:i+2)==' - ' .or. string(i:i+2)==' = ') then
                this%num_species=this%num_species+1
            else
                continue
            end if
        end do
        call this%allocate_reaction()
        j=1
        k=ind_init
        i=ind_init
        do while (j<this%num_species)
            if (string(i:i+2)==' + ' .or. string(i:i+2)==' - ') then
                call this%species(j)%set_name(string(k:i-1))
                !print *, this%species(j)%name
                if (string(i-1:i-1)=='+') then
                    call this%species(j)%set_valence(1)
                else if (string(i-1:i-1)=='-') then
                    call this%species(j)%set_valence(-1)
                else if (string(i-2:i-2)=='+' .or. string(i-2:i-2)=='-') then
                    read(string(i-1:i-1),*) this%species(j)%valence
                    !call this%species(j)%set_valence(1)
                end if
                k=i+3
                j=j+1
            end if
            if (i<n-2) then
                i=i+1
            else
                exit
            end if
        end do
        plus_ind=index(string(eq_ind+1:n),'+')+eq_ind
        minus_ind=index(string(eq_ind+1:n),'-')+eq_ind
        if (plus_ind>eq_ind .and. minus_ind>eq_ind) then
            sign_ind=min(plus_ind,minus_ind)
        else if (plus_ind>eq_ind) then
            sign_ind=plus_ind
        else if (minus_ind>eq_ind) then
            sign_ind=minus_ind
        else
            sign_ind=n+2
        end if
        if (present(defined_species)) then
            call defined_species%set_name(string(eq_ind+2:sign_ind-2))
        end if
        call this%species(j)%set_name(string(eq_ind+2:sign_ind-2))
        j=j+1
        if (sign_ind<n+2) then
            k=sign_ind+2
            i=sign_ind+2
            do while (j<this%num_species)
                if (string(i:i+2)==' + ' .or. string(i:i+2)==' - ') then
                    call this%species(j)%set_name(string(k:i-1))
                    !print *, this%species(j)%name
                    if (string(i-1:i-1)=='+') then
                        call this%species(j)%set_valence(1)
                    else if (string(i-1:i-1)=='-') then
                        call this%species(j)%set_valence(-1)
                    else if (string(i-2:i-2)=='+' .or. string(i-2:i-2)=='-') then
                        read(string(i-1:i-1),*) this%species(j)%valence
                        !call this%species(j)%set_valence(1)
                    end if
                    k=i+3
                    j=j+1
                end if
                if (i<n-2) then
                    i=i+1
                else
                    exit
                end if
            end do
            call this%species(j)%set_name(string(k:n))
            if (string(n-1:n-1)=='+') then
                call this%species(j)%set_valence(1)
            else if (string(n-1:n-1)=='-') then
                call this%species(j)%set_valence(-1)
            else if (string(n-2:n-2)=='+' .or. string(n-2:n-2)=='-') then
                read(string(n-1:n-1),*) this%species(j)%valence
            end if
        end if
    end if
end subroutine