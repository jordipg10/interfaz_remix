function mass_balance_error_ADE_trans_PMF_discharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
    use transport_transient_m, only: transport_1D_transient_c
    implicit none
    class(transport_1D_transient_c), intent(in) :: this
    real(kind=8), intent(in) :: conc_old(:)
    real(kind=8), intent(in) :: conc_new(:)
    real(kind=8), intent(in) :: Delta_t
    real(kind=8), intent(in) :: Delta_x
    real(kind=8) :: mass_bal_err
    
    real(kind=8) :: sum_time,sum_r,inf_term,out_term,dummy_inf,D,q_1
    integer(kind=4) :: i,n
    
    n=size(conc_old)
    
    sum_time=0d0
    
    do i=1,n
        sum_time=sum_time+this%tpt_props_heterog%porosity(i)*(conc_new(i)-conc_old(i))
    end do
    sum_time=sum_time*Delta_x
    
    sum_r=0d0
    do i=1,n
        sum_r=sum_r+this%tpt_props_heterog%source_term(i)*conc_old(i)
    end do
    sum_r=sum_r*Delta_x*Delta_t
    
    inf_term = this%BCs%flux_inf*this%BCs%conc_inf
    out_term = - this%BCs%flux_out*conc_old(n)
    
    inf_term=inf_term*Delta_t
    out_term=out_term*Delta_t
    
    mass_bal_err=sum_time-sum_r-inf_term-out_term
end function