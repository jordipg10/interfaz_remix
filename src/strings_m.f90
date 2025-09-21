module strings_m
    interface
        subroutine compare_str_arrays(str_arr_1,str_arr_2,flag,indices)
            implicit none
            character(len=*), intent(in) :: str_arr_1(:)
            character(len=*), intent(in) :: str_arr_2(:)
            logical, intent(out) :: flag
            integer(kind=4), intent(out), allocatable, optional :: indices(:)
        end subroutine
    end interface

    contains

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
end module