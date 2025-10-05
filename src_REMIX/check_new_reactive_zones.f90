!> This subroutine checks if there are new reactive zones due to changes in concentrations
subroutine check_new_reactive_zones(this,i,tolerance)
    use chemistry_m, only: chemistry_c, reactive_zone_c
    use vectors_m, only: inf_norm_vec_int
    implicit none
    class(chemistry_c) :: this
    integer(kind=4), intent(in)  :: i !> reactive zone index
    real(kind=4), intent(in) :: tolerance !> for concentration of non-flowing species
    
    integer(kind=4) :: j,k,l,p,r,s,m,ii,kk,ll,old_num_react_zones,num_new_react_zones,num_sol_conc_zero,num_nf_sp_new,num_tar_sol,&
    num_rep_rows
    integer(kind=4), allocatable :: tar_sol_indices(:),old_conc_zero_flag(:,:),new_conc_zero_flag(:,:),ind_sol_conc_zero(:),&
    ind_conc_zero(:,:),ind_rep_rows(:)
    type(reactive_zone_c), allocatable :: new_react_zones(:),old_react_zones(:)
    
    old_react_zones=this%reactive_zones !> old reactive zones
    old_num_react_zones=this%num_reactive_zones !> old number of reactive zones
    
    call this%link_target_solids_reactive_zone(i,tar_sol_indices) !> we get the indices of targets solids ('tar_sol_indices') associated to this reactive zone
    num_tar_sol=size(tar_sol_indices) !> number of target solids associated to this reactive zone
    num_sol_conc_zero=0 !> number of target solids with at least one solid with zero concentration
    
!> We count the number of target solids with at least one zero concentration solid
    do j=1,num_tar_sol
        do k=1,this%reactive_zones(i)%num_non_flowing_species
            if (abs(this%target_solids(tar_sol_indices(j))%concentrations(k))<tolerance) then
                num_sol_conc_zero=num_sol_conc_zero+1
                exit
            end if
        end do
    end do
    allocate(ind_sol_conc_zero(num_sol_conc_zero)) !> indices of target solids with at least one zero concentration solid
    allocate(old_conc_zero_flag(num_sol_conc_zero,this%reactive_zones(i)%num_non_flowing_species))
    l=0
    do j=1,num_tar_sol
        do k=1,this%reactive_zones(i)%num_non_flowing_species
            if (abs(this%target_solids(tar_sol_indices(j))%concentrations(k))<tolerance) then
                l=l+1
                ind_sol_conc_zero(l)=tar_sol_indices(j)
                if (l==num_sol_conc_zero) exit
            end if
        end do
    end do
    do l=1,num_sol_conc_zero
        do k=1,this%reactive_zones(i)%num_non_flowing_species
            if (abs(this%target_solids(ind_sol_conc_zero(l))%concentrations(k))<tolerance) then
                old_conc_zero_flag(l,k)=1
            else
                old_conc_zero_flag(l,k)=0
            end if
        end do
    end do
    num_rep_rows=0
    do l=1,num_sol_conc_zero-1
        do m=l+1,num_sol_conc_zero
            if (inf_norm_vec_int(old_conc_zero_flag(l,:)-old_conc_zero_flag(m,:))<tolerance) then
                num_rep_rows=num_rep_rows+1
                exit
            end if
        end do
    end do
    allocate(ind_rep_rows(num_rep_rows))
    p=0
    do l=1,num_sol_conc_zero-1
        do m=l+1,num_sol_conc_zero
            if (inf_norm_vec_int(old_conc_zero_flag(l,:)-old_conc_zero_flag(m,:))<tolerance) then
                p=p+1
                ind_rep_rows(p)=l
                exit
            end if
        end do
    end do
    num_new_react_zones=num_sol_conc_zero-num_rep_rows
    allocate(new_react_zones(num_new_react_zones))
    allocate(new_conc_zero_flag(num_new_react_zones,this%reactive_zones(i)%num_non_flowing_species))
    p=1
    r=0
    do l=1,num_sol_conc_zero
        if (l/=ind_rep_rows(p)) then
            r=r+1
            new_conc_zero_flag(r,:)=old_conc_zero_flag(l,:)
            if (p<num_rep_rows) then
                p=p+1
            else
                exit
            end if
        end if
    end do
    do s=1,num_sol_conc_zero-ind_rep_rows(num_rep_rows)
        r=r+1
        new_conc_zero_flag(r,:)=old_conc_zero_flag(ind_rep_rows(num_rep_rows)+s,:)
    end do
    num_nf_sp_new=0
    do ii=1,num_new_react_zones
        num_nf_sp_new=this%reactive_zones(i)%num_non_flowing_species-sum(new_conc_zero_flag(ii,:))
        call new_react_zones(ii)%allocate_non_flowing_species(num_nf_sp_new)
        call new_react_zones(ii)%set_chem_syst_react_zone(this%chem_syst)
        ll=0
        do kk=1,this%reactive_zones(i)%num_non_flowing_species
            if (new_conc_zero_flag(ii,kk)==0) then
                ll=ll+1
                call new_react_zones(ii)%non_flowing_species(ll)%assign_species(this%reactive_zones(i)%non_flowing_species(kk))
            end if
        end do
    end do
        
    deallocate(this%reactive_zones)
    call this%allocate_reactive_zones(old_num_react_zones+num_new_react_zones)
        
    do l=1,old_num_react_zones
        this%reactive_zones(l)=old_react_zones(l)
    end do
    do l=1,num_new_react_zones
        this%reactive_zones(old_num_react_zones+l)=new_react_zones(l)
    end do
        
    deallocate(tar_sol_indices,old_conc_zero_flag,new_conc_zero_flag,ind_sol_conc_zero,new_react_zones,old_react_zones)    
end subroutine