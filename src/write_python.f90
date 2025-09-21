!> Writes data and results of 1D reactive transport object necessary for Python unit tests
subroutine write_python(this,path)
    use RT_1D_m, only: RT_1D_c, RT_1D_transient_c
    
    implicit none
    class(RT_1D_c), intent(in) :: this                          !> 1D reactive transport object
    character(len=*), intent(in) :: path                        !> path output
    
    integer(kind=4) :: i,j
        
    select type(this)
    type is (RT_1D_transient_c)
        ! open(9988,file=trim(path)//'u_init.dat')
        ! do i=1,this%chemistry%chem_syst%speciation_alg%num_aq_prim_species
        !     write(9988,"(*(ES15.5))") (dot_product(&
        !     this%chemistry%target_waters_init(j)%solid_chemistry%reactive_zone%speciation_alg%comp_mat(i,:),&
        !     this%chemistry%target_waters_init(j)%get_conc_nc()), j=1,this%chemistry%num_target_waters)
        ! end do
        ! close(9988)
        ! open(999,file=trim(path)//'c1_aq_init.dat')
        ! do i=1,this%chemistry%chem_syst%speciation_alg%num_aq_prim_species
        !     write(999,"(*(ES15.5))") (this%chemistry%target_waters_init(j)%concentrations(i), j=1,this%chemistry%num_target_waters)
        ! end do
        ! close(999)
        ! open(9999,file=trim(path)//'gamma_aq_init.dat')
        ! do i=1,this%chemistry%chem_syst%aq_phase%num_species
        !     write(9999,"(*(ES15.5))") (10**(this%chemistry%target_waters_init(j)%log_act_coeffs(i)), &
        !     j=1,this%chemistry%num_target_waters)
        ! end do
        ! close(9999)
        open(998,file=trim(path)//'lambdas_filas.dat')
        if (this%transport%time_discr%int_method.eq.1) then
            write(998,"(20x,*(ES20.10))") this%transport%X_mat%diag(1), this%transport%X_mat%super(1)
            do i=2,this%transport%mixing_ratios_conc%num_cols-1
                write(998,"(*(ES20.10))")  this%transport%X_mat%sub(i-1), this%transport%X_mat%diag(i), &
                    this%transport%X_mat%super(i)
            end do
            write(998,"(*(ES20.10))") this%transport%X_mat%sub(this%transport%mixing_ratios_conc%num_cols-1), &
                this%transport%X_mat%diag(this%transport%mixing_ratios_conc%num_cols)
        else if (this%transport%time_discr%int_method.eq.2) then
            do i=1,this%transport%mixing_ratios_conc%num_cols
                write(998,"(*(ES20.10))") (this%transport%mixing_ratios_mat_conc_dom(j,i), &
                    j=1,this%transport%mixing_ratios_conc%num_cols)
            end do
        end if
        close(998)
        ! open(99855,file=trim(path)//'lambdas_Rk_filas.dat')
        ! if (this%transport%time_discr%int_method.eq.1) then
        !     write(998,"(20x,*(ES20.10))") this%transport%X_mat%diag(1), this%transport%X_mat%super(1)
        !     do i=2,this%transport%mixing_ratios_conc%num_cols-1
        !         write(99855,"(*(ES20.10))")  this%transport%X_mat%sub(i-1), this%transport%X_mat%diag(i), &
        !             this%transport%X_mat%super(i)
        !     end do
        !     write(99855,"(*(ES20.10))") this%transport%X_mat%sub(this%transport%mixing_ratios_conc%num_cols-1), &
        !         this%transport%X_mat%diag(this%transport%mixing_ratios_conc%num_cols)
        ! else if (this%transport%time_discr%int_method.eq.2) then
        !     do i=1,this%transport%mixing_ratios_Rk_init%num_cols
        !         write(99855,"(*(ES20.10))") (this%transport%A_mat_inv(j,i), &
        !             j=1,this%transport%mixing_ratios_Rk_init%num_cols)
        !     end do
        ! end if
        ! close(99855)
        ! open(997,file=trim(path)//'lambdas_cols.dat')
        ! if (this%transport%time_discr%int_method.eq.1) then
        !     do i=1,this%transport%mixing_ratios_conc%num_cols
        !         write(997,"(*(ES20.10))") (this%transport%mixing_ratios_conc%cols(i)%col_1(j), &
        !             j=1,this%transport%mixing_ratios_conc%cols(i)%dim)
        !     end do
        ! else if (this%transport%time_discr%int_method.eq.2) then
        !     do i=1,this%transport%mixing_ratios_conc%num_cols
        !         write(997,"(*(ES20.10))") (this%transport%mixing_ratios_mat_conc_dom(i,j), &
        !             j=1,this%transport%mixing_ratios_conc%num_cols)
        !     end do
        ! end if
        ! close(997)
        ! open(9977,file=trim(path)//'rk_init.dat')
        ! if (this%chemistry%chem_syst%num_kin_reacts>0) then
        !     do i=1,this%chemistry%chem_syst%num_aq_kin_reacts
        !         write(9977,"(*(ES20.10))") (this%chemistry%target_waters_init(j)%rk(i), &
        !             j=1,this%chemistry%num_target_waters)
        !     end do
        ! end if
        ! close(9977)
        ! open(99777,file=trim(path)//'rk.dat')
        ! if (this%chemistry%chem_syst%num_kin_reacts>0) then
        !     do i=1,this%chemistry%chem_syst%num_aq_kin_reacts
        !         write(99777,"(*(ES20.10))") (this%chemistry%target_waters(j)%rk(i), &
        !             j=1,this%chemistry%num_target_waters)
        !     end do
        ! end if
        ! close(99777)
        ! open(12,file=trim(path)//'conc_nc_aq_init.dat')
        ! !if (this%chemistry%chem_syst%num_kin_reacts>0) then
        !     do i=1,this%chemistry%chem_syst%speciation_alg%num_aq_var_act_species
        !         write(12,"(*(ES20.10))") (this%chemistry%target_waters_init(j)%concentrations(&
        !             this%chemistry%target_waters_init(j)%ind_var_act_species(i)), &
        !             j=1,this%chemistry%num_target_waters)
        !     end do
        ! !end if
        ! close(12)
        ! open(11,file=trim(path)//'conc_nc_aq.dat')
        ! !if (this%chemistry%chem_syst%num_kin_reacts>0) then
        !     do i=1,this%chemistry%chem_syst%speciation_alg%num_aq_var_act_species
        !         write(11,"(*(ES20.10))") (this%chemistry%target_waters(j)%concentrations(&
        !             this%chemistry%target_waters_init(j)%ind_var_act_species(i)), &
        !             j=1,this%chemistry%num_target_waters)
        !     end do
        ! !end if
        ! close(11)
    end select
end subroutine write_python
!> End of file write_python.f90