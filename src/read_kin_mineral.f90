!> Reads parameters for mineral kinetic reaction
subroutine read_kin_mineral(this)
    use kin_mineral_m
    use strings_m
    implicit none
    class(kin_mineral_c) :: this
    
    character(len=256) :: str
            
    open(unit=1,file='kinetics.dat',status='old',action='read')
    do
        read(1,*) str
        if (str=='end') exit
        if (str==this%species(1)%name) then
            backspace(1)
            exit
        else
            continue
        end if
    end do
    close(1)
end subroutine