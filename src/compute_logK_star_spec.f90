subroutine compute_logK_star(this,K)
    use speciation_algebra_m, only: speciation_algebra_c
    implicit none
    class(speciation_algebra_c) :: this
    real(kind=8), intent(in) :: K(:) ! eq csts
   
    this%logK_star=matmul(this%inv_Se_nc_2,log10(K))
end subroutine