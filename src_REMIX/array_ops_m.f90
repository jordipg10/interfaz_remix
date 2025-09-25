!> Array operations module
module array_ops_m
    implicit none
    save
    contains
        subroutine append_int_1D_array(array,new_elem) !> appends element in integer vector
            implicit none
            integer(kind=4), intent(inout), allocatable :: array(:)
            integer(kind=4), intent(in) :: new_elem
            
            integer(kind=4) :: i
            integer(kind=4), allocatable :: aux_array(:)
            
            !print *, size(array)
            aux_array=array
            if (allocated(array)) then
                deallocate(array)
            end if
            allocate(array(size(aux_array)+1))
            do i=1,size(array)-1
                array(i)=aux_array(i)
            end do
            array(size(array))=new_elem
        end subroutine
        
        subroutine is_int_in_1D_array(int,array,flag,ind)
        !> checks if integer is in array
            implicit none
            integer(kind=4), intent(in) :: int !> integer to check
            integer(kind=4), intent(in) :: array(:) !> array of integers to check
            logical, intent(out) :: flag !> TRUE if integer is in array, FALSE otherwise
            integer(kind=4), intent(out), optional :: ind !> index of integer in array
            
            integer(kind=4) :: i
            
            flag=.false.
            if (present(ind)) then
                ind=0
            end if
            do i=1,size(array)
                if (array(i)==int) then
                    flag=.true.
                    if (present(ind)) then
                        ind=i
                    end if
                    exit
                else
                    continue
                end if
            end do
        end subroutine
        
        function get_ind_int(int,array) result(ind)
            integer(kind=4), intent(in) :: int
            integer(kind=4), intent(in) :: array(:)
            integer(kind=4) :: ind
            
            integer(kind=4) :: i
            do i=1,size(array)
                if (array(i)==int) then
                    ind=i
                end if
            end do
        end function
        
end module