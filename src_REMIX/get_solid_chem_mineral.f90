!subroutine get_solid_chem_mineral(this,aq_species,activities,act_cat,aq_species_ind)!,react_surf,temp)
!    use kin_mineral_m
!    implicit none
!    class(kin_mineral_c), intent(in) :: this
!    class(aq_species_c), intent(in) :: aq_species(:)
!    real(kind=8), intent(in) :: activities(:) !> aqueous species activities
!    real(kind=8), intent(out) :: act_cat(:) !> activities catalysers
!    integer(kind=4), intent(out), optional :: aq_species_ind(:) !> indices of aqueous species of reaction
!    !real(kind=8), intent(out) :: react_surf
!    !real(kind=8), intent(out) :: temp
!            
!    integer(kind=4) :: i,j,k,l,m,p,n,DOC_ind
!    real(kind=8), parameter :: epsilon=1d-6
!    type(species_c) :: DOC !> CH2O
!    logical :: DOC_flag
!    
!    i=1
!    j=1
!    !k=1
!    do
!        if (aq_species(i)%name==this%species(j)%name) then
!            aq_species_ind(j)=i
!            !if (this%species(j)%name==this%params%catalysers(k)%name) then
!            !>    act_cat(k)=activities(i)
!            !>    if (k<this%params%num_cat) then
!            !>        k=k+1
!            !>    end if
!            !end if
!            if (j<this%num_species) then
!                j=j+1
!                i=1
!            else
!                exit
!            end if
!        else if (i<size(aq_species)) then
!            i=i+1
!        else
!            exit
!        end if
!    end do
!    
!    i=1
!    k=1
!    do
!        if (aq_species(i)%name==this%params%catalysers(k)%name) then
!            act_cat(k)=activities(i)
!            !if (this%species(j)%name==this%params%catalysers(k)%name) then
!            !>    act_cat(k)=activities(i)
!            !>    if (k<this%params%num_cat) then
!            !>        k=k+1
!            !>    end if
!            !end if
!            if (k<this%params%num_cat) then
!                k=k+1
!                i=1
!            else
!                exit
!            end if
!        else if (i<size(aq_species)) then
!            i=i+1
!        else
!            exit
!        end if
!    end do
!            
!   !> l=1
!   !> m=1
!   !> n=1
!   !> if (size(conc)/=size(species)) error stop "Dimension error in get_conc_kin_mineral"
!   !>                 
!   !error stop "Subroutine 'get_conc_kin_mineral' not implemented yet"
!end subroutine