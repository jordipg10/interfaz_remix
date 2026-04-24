!> \file array_ops_m.f90
!> \brief Array operations utility module for integer and real array manipulations.
!> \details
!> Provides helper subroutines and functions for common array operations used
!> throughout the RT_code codebase, including:
!> - Appending elements to dynamically-allocated integer arrays
!> - Searching for integers within arrays (membership tests)
!> - Retrieving the index of an integer in an array
!>
!> These routines are used extensively for managing index arrays, species lists,
!> and other dynamic data structures in the reactive transport framework.
!>
!> \author Jordi
!> \date Unknown
!> \ingroup algebra

!> \brief Array operations utility module.
!> \details Contains general-purpose array manipulation subroutines and functions.
module array_ops_m
    implicit none
    save
    contains

        !> \brief Append an integer element to a 1D allocatable integer array.
        !> \details Resizes the array by one and places the new element at the end.
        !> Uses a temporary copy to preserve existing data during reallocation.
        !> \param[in,out] array  Allocatable integer array to be extended
        !> \param[in]     new_elem  Integer value to append
        subroutine append_int_1D_array(array,new_elem)
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
        
        !> \brief Check if an integer value exists in a 1D integer array.
        !> \details Performs a linear search through the array. Optionally returns
        !> the index of the first match.
        !> \param[in]  int    Integer value to search for
        !> \param[in]  array  Array of integers to search in
        !> \param[out] flag   .true. if the integer is found, .false. otherwise
        !> \param[out] ind    (Optional) Index of the first occurrence in the array (0 if not found)
        subroutine is_int_in_1D_array(int,array,flag,ind)
            implicit none
            integer(kind=4), intent(in) :: int !< Integer value to search for
            integer(kind=4), intent(in) :: array(:) !< Array of integers to search in
            logical, intent(out) :: flag !< .true. if integer is found, .false. otherwise
            integer(kind=4), intent(out), optional :: ind !< Index of integer in array (0 if not found)
            
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
        
        !> \brief Get the index of an integer in a 1D array.
        !> \details Performs a linear search and returns the index of the matching element.
        !> \warning Undefined behavior if the integer is not present in the array.
        !> \param[in] int    Integer value to search for
        !> \param[in] array  Array to search in
        !> \return Index of the matching element in the array
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