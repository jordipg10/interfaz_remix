!> This module contains the parameters to compute conventional specific volume of solute species
module params_spec_vol_m
    implicit none
    save
    type, public :: params_spec_vol_Redlich_c !> Redlich-type equation (PHREEQC)
        real(kind=8) :: a1 !>  
        real(kind=8) :: a2 !> 
        real(kind=8) :: a3 !> 
        real(kind=8) :: a4 !
        real(kind=8) :: W !> 
        real(kind=8) :: ion_size_param !> Debye-Huckel
        real(kind=8) :: i1 !> 
        real(kind=8) :: i2 !
        real(kind=8) :: i3 !
        real(kind=8) :: i4 !> 
    contains
        procedure, public :: set_params
    end type
    
    contains
        subroutine set_params(this,coeffs)
            implicit none
            class(params_spec_vol_Redlich_c) :: this
            real(kind=8), intent(in) :: coeffs(:) !> must be in same order as PHREEQC
            
            if (size(coeffs)/=10) then
                error stop "Dimension error in coefficients specific volume Redlich"
            else
                this%a1=coeffs(1)
                this%a2=coeffs(2)
                this%a3=coeffs(3)
                this%a4=coeffs(4)
                this%W=coeffs(5)
                this%ion_size_param=coeffs(6)
                this%i1=coeffs(7)
                this%i2=coeffs(8)
                this%i3=coeffs(9)
                this%i4=coeffs(10)
            end if
        end subroutine
end module