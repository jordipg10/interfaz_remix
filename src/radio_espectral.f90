!> Spectral radius
function radio_espectral(lambda) result(rho)
    implicit none
    real(kind=8), intent(in) :: lambda(:)
    real(kind=8) :: rho
    rho=max(maxval(lambda),maxval(-lambda))
end function 