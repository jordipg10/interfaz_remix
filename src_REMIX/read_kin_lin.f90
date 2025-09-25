!> Reads parameters for linear kinetic reaction
subroutine read_kin_lin(this)
    use lin_kin_reaction_m
    use strings_m
    implicit none
    class(lin_kin_reaction_c) :: this
    
    character(len=256) :: str
            
    open(unit=1,file='reacciones_lineales.dat',status='old',action='read')
    do
        read(1,*) str
        if (str=='end') exit
        if (str==this%species(1)%name) then
            backspace(1)
            read(1,*) str, this%lambda
            exit
        else
            continue
        end if
    end do
    close(1)
end subroutine