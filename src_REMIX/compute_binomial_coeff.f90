!> This functions computes binomial coefficient
function compute_binomial_coeff(m,n) result(bin_coeff)
    use special_fcts_m, only: fact
    implicit none
    integer(kind=4), intent(in) :: m !> up
    integer(kind=4), intent(in) :: n !> down
    integer(kind=4) :: bin_coeff
    
    integer(kind=4) :: i,aux,small_fact
    
    if (n>m) error stop "Error in binomial coefficient"
    
    bin_coeff=fact(m)
    bin_coeff=bin_coeff/fact(n)
    bin_coeff=bin_coeff/fact(m-n)
end function