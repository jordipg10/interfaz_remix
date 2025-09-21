!> This subroutine reads the time discretisation used for the WMA from a file ended in "_WMA_discr.dat"
subroutine read_time_discretisation(this,root)
    use RT_1D_m, only: RT_1D_transient_c, RT_1D_c
    use time_discr_m, only: time_discr_homog_c, time_discr_heterog_c, time_discr_c
    use vectors_m, only: inf_norm_vec_real
    use, intrinsic :: iso_fortran_env
    implicit none
    class(RT_1D_c) :: this                                      !> transient 1D reactive transport object
    !integer(kind=4), intent(in) :: unit                         !> unit of file
    character(len=*), intent(in) :: root                        !> root of file
    
    integer(KIND=4) :: Num_time                                 !< number of time steps
    integer(KIND=4) :: int_method_chem_reacts                   !< integration method chemical reactions
    integer(KIND=4) :: unit                                     !< unit of file
    character(len=256) :: label
    type(time_discr_homog_c), target :: time_discr_homog
    type(time_discr_heterog_c), target :: time_discr_heterog
    class(time_discr_c), pointer :: time_discr=>NULL()
    real(kind=8), allocatable :: Delta_t(:)                     !< time step array
    real(kind=8), parameter :: EPSILON = 1d-16        !< machine epsilon for kind=8
    
    unit=5 !> arbitrary unit
    select type (this)
    type is (RT_1D_transient_c)
    !> We open file with time discretisation data
        open(unit,file=root//'_WMA_discr.dat',status='old',action='read')
        do 
            read(unit,*) label
            if (label.eq.'end') then
                exit
            else if (label.eq.'TIME') then
            !> we read number of time steps
                read(unit,*) Num_time
                allocate(Delta_t(Num_time))
            !< we read time step or steps
                read(unit,*) Delta_t(1) !> autentica chapuza
                ! if (inf_norm_vec_real(Delta_t-Delta_t(1))<EPSILON) then !< if time step is uniform, then we use homogeneous subclass
                !     !call time_discr_homog%set_Delta_t_homog(Delta_t(1))
                !     !allocate(time_discr)
                !     !time_discr=>time_discr_homog
                !     allocate(time_discr_homog_c :: time_discr) !> chapuza
                ! else
                !     !call time_discr_heterog%set_Delta_t_heterog(Delta_t) !< if time step is not uniform, then we use heterogeneous subclass
                !     !time_discr=>time_discr_heterog
                !     allocate(time_discr_heterog_c :: time_discr) !> chapuza
                ! end if
                allocate(time_discr_homog_c :: time_discr) !> chapuza
                Delta_t(2:Num_time)=Delta_t(1) !> chapuza
                call time_discr%set_Delta_t(Delta_t)
            !< we read time integration method for chemical reactions
                read(unit,*) int_method_chem_reacts
                call this%set_int_method_chem_reacts(int_method_chem_reacts)
            !< we set the number of time steps attribute
                call time_discr%set_Num_time(Num_time)
            !< we compute the final time
                call time_discr%compute_Final_time()
                call this%transport%set_time_discr(time_discr)
            else
                continue
            end if
        end do
        close(unit)
    end select
end subroutine