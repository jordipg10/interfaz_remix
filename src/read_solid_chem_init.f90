!> subroutine read_solid_chem_init(this,filename,reactive_zones,line,num_tar)
!>    use solid_chemistry_m
!>    implicit none
!>    class(solid_chemistry_c) :: this
!>    character(len=*), intent(in) :: filename
!>    class(reactive_zone_c), intent(in) :: reactive_zones(:) !> from the chemistry class
!>    integer(kind=4), intent(inout) :: line
!>    integer(kind=4), intent(out) :: num_tar !> number of targets associated 
!>    
!>    integer(kind=4) :: j,solid_chem_ind,react_zone_ind,nf_sp_ind
!>    type(species_c) :: species,constrain
!>    character(len=256) :: str
!>    real(kind=8) :: aux,conc,temp
!>    logical :: flag
!>  
!>    read(5,*) react_zone_ind, num_tar, this%temp
!>    if (react_zone_ind<1 .or. react_zone_ind>size(reactive_zones)) error stop
!>    call this%set_reactive_zone(reactive_zones(react_zone_ind))
!>    call this%allocate_concentrations()
!>    j=1
!>    do
!>        line=line+1
!>        read(5,*) species%name
!>        if (species%name=='*') exit
!>        call this%reactive_zone%is_nf_species_in_react_zone(species,flag,nf_sp_ind)
!>        if (flag==.true.) then
!>            backspace(5)
!>            read(5,*) species%name, this%concentrations(j)
!>            if (j<this%reactive_zone%num_non_flowing_species) then
!>                j=j+1
!>            end if
!>        else if (species%name=='') then
!>            continue
!>        else
!>            error stop "This non-flowing species does not belong to the reactive zone"
!>        end if
!>    end do
!>                  
!end subroutine