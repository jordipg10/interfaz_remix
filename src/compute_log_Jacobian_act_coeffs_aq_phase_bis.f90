!!> Computes Jacobian log_10 activity coefficients aqueous species with respect to log_10 concentrations
!subroutine compute_log_Jacobian_act_coeffs_aq_phase(this,d_log_gamma_d_I,d_I_d_c_aq,conc,log_Jacobian_act_coeffs)
!    use aq_phase_m
!    implicit none
!    class(aq_phase_c) :: this
!    real(kind=8), intent(in) :: d_log_gamma_d_I(:) !> derivative of log_10(activity coefficients) with respect to ionic activity
!    real(kind=8), intent(in) :: d_I_d_c_aq(:) !> Jacobian of ionic activity with respect to aqueous cocnentrationa
!    real(kind=8), intent(in) :: conc(:) !> concentration of aqueous species
!    real(kind=8), intent(out) :: log_Jacobian_act_coeffs(:,:) !> (must be allocated)   
!
!    integer(kind=4) :: i,j
!        
!    do i=1,size(d_log_gamma_d_I)
!        do j=1,size(d_log_gamma_d_I)
!            log_Jacobian_act_coeffs(i,j)=d_log_gamma_d_I(i)*d_I_d_c_aq(j)*conc(j)*log(1d1)
!            !> Debye-Huckel:
!                !log_Jacobian_act_coeffs(i,j)=(25d-2*log_act_coeffs(i)*conc(j)*log(1d1)*this%aq_species(j)%valence**2)/(ionic_act)
!            !> Davies:
!                !log_Jacobian_act_coeffs(i,j)=(25d-2*log_act_coeffs(i)*conc(j)*log(1d1)*this%aq_species(j)%valence**2)/(ionic_act*(1d0+sqrt(ionic_act)))
!            !> Debye-Huckel ampliado:
!                !log_Jacobian_act_coeffs(i,j)=(25d-2*log_act_coeffs(i)*conc(j)*log(1d1)*this%aq_species(j)%valence**2)/(ionic_act*(1d0+this%params_act_coeff%ion_size_param*this%params_act_coeff%B*sqrt(ionic_act)))
!        end do
!    end do
!end subroutine