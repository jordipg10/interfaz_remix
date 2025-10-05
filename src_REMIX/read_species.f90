subroutine read_species(this,str)
    use species_m, only: species_c
    implicit none
    class(species_c) :: this
    character(len=*), intent(in) :: str !> we assumed it is already trimmed
    
    integer(kind=4) :: len_str,index_pos,index_neg,valence
    
    len_str=len(str)
    call this%set_name(str)
    index_pos=index(str,'+')
    index_neg=index(str,'-')
    if (index_pos/=0) then
        if (index_pos==len_str-1) then
            valence=1
        else
            read(str(index_pos:len_str-1),*) this%valence
        end if
    else if (index_neg/=0) then
        if (index_neg==len_str-1) then
            valence=-1
        else
            read(str(index_neg:len_str-1),*) this%valence
        end if
    else
        call this%set_valence(0)
    end if
end subroutine