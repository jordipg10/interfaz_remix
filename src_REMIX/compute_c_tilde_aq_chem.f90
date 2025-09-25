!> This subroutine computes the concentration of aqueous variable activity species after mixing
function compute_c_tilde(this,mixing_ratios,conc_old) result(c_tilde)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Argumnets
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(in) :: mixing_ratios(:)  !> the first element corresponds to argument "this", and the rest correspond to argument "mixing waters" (in the same order)
    real(kind=8), intent(in) :: conc_old(:,:) !> concentration of aqueous variable activity species before mixing
    real(kind=8), allocatable :: c_tilde(:) !> concentration of aqueous variable activity species after mixing
!> Variables
    integer(kind=4) :: i !> index aqueous variable activity species
    integer(kind=4) :: j !> index mixing waters
!> Pre-process
    allocate(c_tilde(size(conc_old,1)))
!> Process
    do i=1,size(conc_old,1)
        c_tilde(i)=mixing_ratios(1)*conc_old(i,1)
        do j=1,size(conc_old,2)-1
            c_tilde(i)=c_tilde(i)+mixing_ratios(j+1)*conc_old(i,j+1)
        end do
    end do
end function 