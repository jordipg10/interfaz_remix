!> Special functions module
module special_fcts_m
    implicit none
    save
    interface
        ! function compute_binomial_coeff(m,n) result(bin_coeff)
        !     implicit none
        !     integer(kind=4), intent(in) :: m !> up
        !     integer(kind=4), intent(in) :: n !> down
        !     integer(kind=4) :: bin_coeff
        ! end function
    end interface
    
    contains
        recursive function fact(n) result(fact_result) !> factorial
            implicit none
            integer(kind=4), intent(in) :: n
            integer(kind=4) :: fact_result
            if (n<0) error stop "Factorial must be non-negative"
            if (n.eq.0) then
                fact_result=1
            else 
                fact_result=n*fact(n-1)
            end if
        end function
        
        function Pochhammer(x,n)
            implicit none
            real(kind=8), intent(in) :: x
            integer(kind=4), intent(in) :: n
            real(kind=8) :: Pochhammer
            Pochhammer=gamma(x+n)/gamma(x)
        end function
        
        function confl_hypergeo_fct_1st_kind(a,b,z) result(F1)
            implicit none
            real(kind=8), intent(in) :: a
            real(kind=8), intent(in) :: b
            complex(kind=8), intent(in) :: z
            real(kind=8) :: F1
            
            integer(kind=4) :: i,n
            
            n=1000000
            F1=0d0
            do i=1,n
                F1=F1+Pochhammer(a,i)*(z**i)/(Pochhammer(b,i)*fact(i))
            end do
        end function
        
        function parab_cyl_fct(nu,z) result(D)
            implicit none
            real(kind=8), intent(in) :: nu
            complex(kind=8), intent(in) :: z
            real(kind=8) :: D
            
            complex(kind=8) :: i=(0d0,1d0)
            
            D=(2**(nu/2)*exp((-z**2)/4)*((-i*z)**(0.25))*(i*z)**(0.25))*confl_hypergeo_fct_1st_kind(-0.5*nu,5d-1,(z**2)/2)/sqrt(z)
        end function

        !> This functions computes binomial coefficient
function compute_binomial_coeff(m,n) result(bin_coeff)
  
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


end module