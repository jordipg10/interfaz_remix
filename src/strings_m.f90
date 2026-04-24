!> \file strings_m.f90
!> \brief String utility module for character array operations.
!> \details
!> Provides helper routines for searching and comparing character string arrays.
!> Used throughout the codebase for name-based lookups of species, reactions, etc.
!>
!> \author Jordi
!> \date Unknown
!> \ingroup algebra
module strings_m
    interface
        !> \brief Compare two string arrays for element-wise equality.
        !> \param[in]  str_arr_1 First string array
        !> \param[in]  str_arr_2 Second string array
        !> \param[out] flag      .true. if arrays match element-wise
        !> \param[out] indices   (Optional) Mapping indices
        subroutine compare_str_arrays(str_arr_1,str_arr_2,flag,indices)
            implicit none
            character(len=*), intent(in) :: str_arr_1(:)
            character(len=*), intent(in) :: str_arr_2(:)
            logical, intent(out) :: flag
            integer(kind=4), intent(out), allocatable, optional :: indices(:)
        end subroutine
    end interface

    contains

    !> \brief Check if a string belongs to a string array.
    !> \details Linear search returning the first matching index.
    !> \param[in]  string  String to search for
    !> \param[in]  array   Array of strings to search in
    !> \param[out] flag    .true. if found
    !> \param[out] index   (Optional) Index of the first match
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