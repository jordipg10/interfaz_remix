!> Solves reactive transport problem. Computes state variables and reaction rates.
subroutine solve_RT_1D(this,root,unit)
    use RT_1D_m, only: RT_1D_c, RT_1D_transient_c
    implicit none
    class(RT_1D_c) :: this                                      !> 1D reactive trans
    character(len=*), intent(in) :: root                        !> root of file with output data
    integer(kind=4), intent(in) :: unit                         !> unit of file with output data
    
    select type (this)
    class is (RT_1D_transient_c)
        !> falta el solver de transporte
        !call this%chemistry%solve_reactive_mixing(root,unit,this%transport%mixing_ratios,this%transport%mixing_waters_indices,&
        !    this%transport%F_mat%diag,this%transport%time_discr,this%int_method_chem_reacts)
    end select
end subroutine solve_RT_1D