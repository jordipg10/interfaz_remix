!program main
!>    !use BCs_subroutines_m
!>    use transport_transient_m
!>    implicit none
!
!>    !> Variables
!>    class(PDE_1D_c), pointer :: my_model=>null()
!>    type(diffusion_1D_transient_c), target :: my_diff_transient
!>    type(transport_1D_transient_c), target :: my_tpt_transient
!>    
!>    !real(kind=8), external :: dot
!>    
!>    integer(kind=4) :: eqn_flag,opcion,info
!>    real(kind=8) :: res
!!****************************************************************************************************************************************************
!>    eqn_flag=2
!>    if (eqn_flag.eqv.1) then
!>        my_model=>my_diff_transient
!>    else
!>        my_model=>my_tpt_transient
!>        !call my_tpt_transient%main_PDE()
!>    end if
!>    opcion=3
!>    call my_model%main_PDE(opcion)
!!!> Lapack
!!>    !call sterf(my_model%trans_mat%diag,my_model%trans_mat%sub,info)
!!>    res = dot(my_model%trans_mat%diag,my_model%trans_mat%diag)
!!>    print *, res
!end program