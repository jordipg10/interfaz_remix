subroutine str_belongs_to(string,array,flag,index)
    implicit none
    character(len=*), intent(in) :: string
    character(len=*), intent(in) :: array(:)
    logical, intent(out) :: flag
    integer(kind=4), intent(out), optional :: index
    
    integer(kind=4) :: i,n
    
    flag=.false.
    n=size(array)
    if (n>0) then
        i=1
        do
            if (string.eq.array(i)) then
                flag=.true.
                index=i
                exit
            else if (i<n) then
                i=i+1
            else
                exit
            end if
        end do
    end if
end subroutine