!subroutine set_stab_params_RT(this,stab_params_obj)
!>    use RT_1D_m
!>    use PDE_transient_m
!>    use transport_stab_params_homog_m
!>    implicit none
!>    class(RT_1D_c) :: this
!>    class(stab_params_c), intent(in), target, optional :: stab_params_obj
!>    select type (this)
!>    class is (RT_1D_transient_c)
!>        if (present(stab_params_obj)) then
!>            this%stab_params=>stab_params_obj
!>        else
!>            select type (tpt=>this%transport)
!>            class is (PDE_1D_transient_c)
!>                this%stab_params=>tpt%stab_params
!>                select type (stab_params=>tpt%stab_params)
!>                type is (stab_params_tpt_homog_c)
!>                    print *, stab_params%alpha, stab_params%Peclet, stab_params%beta
!>                end select
!>            end select
!>        end if
!>    end select
!>      
!end subroutine