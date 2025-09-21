subroutine eliminate_non_flowing_species(this,nf_ind)
    use reactive_zone_Lagr_m
    implicit none
    class(reactive_zone_c) :: this
    integer(kind=4), intent(in) :: nf_ind(:) ! indices of non-flowing species to be eliminated from reactive zone
    
    integer(kind=4) :: i,j,k,l,flag,new_num_non_flowing_species
    integer(kind=4), allocatable :: solid_ind(:),old_min_ind(:)
    type(reactive_zone_c) :: new_react_zone
    type(species_c), allocatable :: new_non_flowing_species(:)

    new_num_non_flowing_species=this%num_non_flowing_species-size(nf_ind)
    allocate(new_non_flowing_species(new_num_non_flowing_species))
    l=1
    j=1
    do
        k=1
        do
            if (k/=nf_ind(j)) then
                !new_react_zone%non_flowing_species(l)=this%non_flowing_species(k)
                !l=l+1
                !if (l<new_react_zone%num_non_flowing_species) then
                !    l=l+1
                !end if
                if (j<size(nf_ind)) then
                    j=j+1
                else
                    new_non_flowing_species(l)=this%non_flowing_species(k)
                    if (l<new_num_non_flowing_species) then
                        l=l+1
                    else
                        exit
                    end if
                end if
            else if (k<this%num_non_flowing_species) then
                k=k+1
            else
                exit
            end if
        end do
    end do
            
         
end subroutine