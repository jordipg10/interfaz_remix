!> Computes factorial of integer
recursive function fact(n) result(res)
    implicit none
    integer(kind=4), intent(in) :: n
    integer(kind=4) :: res
    
    if (n<0) error stop "Factorial must be non-negative"
    if (n == 0) then
        res=1
    else 
        res=n*fact(n-1)
    end if
end function