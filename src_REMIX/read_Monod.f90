!!> Reads Monod redox reaction from database 'reacciones_monod.dat' (created by me)
!subroutine read_Monod(this)
!    use Monod_m
!    implicit none
!    class(redox_kin_c) :: this
!            
!    integer(kind=4) :: i,j,n_t,index,water_flag
!    character(len=256) :: str
!    
!    open(unit=1,file='reacciones_monod.dat',status='old',action='read')
!    do
!        read(1,*) str
!        !print *, str
!        if (str=='end') exit
!        if (str==this%name) then
!            read(1,*) this%params%n_inh
!            if (this%params%n_inh>0) then
!                call this%params%allocate_inhibitors()
!                backspace(1)
!                read(1,*) this%params%n_inh, (this%params%inhibitors(j)%name, this%params%k_inh(j), j=1,this%params%n_inh)
!            end if
!            read(1,*) this%params%n_M
!            call this%params%allocate_TEAs()
!            backspace(1)
!            read(1,*) this%params%n_M, (this%params%TEAs(j)%name, this%params%k_M(j), j=1,this%params%n_M)
!            call this%params%compute_n_t()
!            read(1,*) this%num_species
!            call this%allocate_reaction()
!            backspace(1)
!            read(1,*) this%num_species, (this%species(j)%name, this%stoichiometry(j), j=1,this%num_species), this%params%rate_cst
!            exit
!        else
!            continue
!        end if
!    end do
!    close(1)
!        
!end subroutine