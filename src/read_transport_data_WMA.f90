!> This subroutine reads the mixing ratios and the mixing water indices from the file '_WMA_lambdas.dat'
!! We assume that the mixing waters are already ordered in upstream & downstream, respectively
!! It is possible to deduce the indices of upstream & downstream target waters from the mixing ratios. This option is not implemented yet.
subroutine read_transport_data_WMA(this,root)
    use BCs_m, only: BCs_t
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c 
    use time_discr_m, only: time_discr_homog_c
    use transport_transient_m, only: transport_1D_transient_c
    use transport_properties_heterog_m, only: tpt_props_heterog_c
    use matrices_m, only: real_array_c, diag_matrix_c

    class(transport_1D_transient_c) :: this !> transport object
    !character(len=*), intent(in) :: path
    !integer(kind=4), intent(in) :: unit !> unit number
    character(len=*), intent(in) :: root !> root name of the file
    !character(len=*), intent(in) :: file_tpt
    
    character(len=256) :: str,label
    character(len=:), allocatable :: str_trim
    type(mesh_1D_Euler_homog_c), target :: mesh
    type(time_discr_homog_c) :: time_discr
    type(BCs_t) :: BCs
    type(tpt_props_heterog_c) :: tpt_props
    integer(kind=4) :: i,j,BCs_label(2),comm_ind,pos,tar_dim,mix_wat_ind,dim,unit
    integer(kind=4), allocatable :: num_mix_waters(:),mix_wat_indices(:)
    real(kind=8) :: Delta_x,Delta_t,r,phi,D
    type(real_array_c) :: mixing_ratios_conc
    logical :: flag_props
    type(diag_matrix_c) :: F_mat
    
    !type(transport_1D_transient_c) :: this
    unit=5 !> arbitrary unit number
    
    
    open(unit,file=root//'_WMA_lambdas.dat',status='old',action='read')
    do 
        read(unit,*) label
        if (label.eq.'end') then
            rewind(unit)
            exit
        else if (label.eq.'MIXING RATIOS') then
            i=0 !> counter target waters
            do
                read(unit,*) tar_dim !> number of target waters
                if (tar_dim.eq.0) then
                    exit
                else
                    i=i+1
                end if
            end do
            allocate(mesh_1D_Euler_homog_c :: this%spatial_discr) !> chapuza
            !this%spatial_discr => mesh !> chapuza
            call this%spatial_discr%set_Num_targets(i)
            call mixing_ratios_conc%allocate_array(i)
            call this%mixing_ratios_Rk%allocate_array(i)
            call this%mixing_ratios_Rk_init%allocate_array(i)
        else
            continue
        end if
    end do
    do 
        read(unit,*) label
        if (label.eq.'end') then
            rewind(unit)
            exit
        ! else if (label.eq.'TRANSPORT PROPERTIES') then
        !     tpt_props%homog_flag=.true. !> homogeneous transport properties (by default)
        !     read(unit,*) flag_props, phi !> flag for homogeneous transport properties and porosity
        !     if (flag_props.eqv..true.) then
        !         allocate(tpt_props%porosity(this%spatial_discr%Num_targets),this%F_mat%diag(this%spatial_discr%Num_targets))
        !         tpt_props%porosity=phi !> we assume that the porosity is the same for all the target waters
        !         this%F_mat%diag=phi  !> storage matrix
        !     else
        !         tpt_props%homog_flag=.false. !> heterogeneous transport properties
        !     end if
        !     call this%set_tpt_props_heterog_obj(tpt_props) !> we set the transport properties object
        else if (label.eq.'MIXING RATIOS') then
            i=0 !> counter targets
            do
                i=i+1
                if (i<=this%spatial_discr%Num_targets) then
                    read(unit,*) dim !> number of mixing ratios in this target
                    call mixing_ratios_conc%cols(i)%allocate_vector(dim)
                    call this%mixing_ratios_Rk%cols(i)%allocate_vector(dim)
                    call this%mixing_ratios_Rk_init%cols(i)%allocate_vector(dim)
                else
                    exit
                end if
            end do
        else
            continue
        end if
    end do
    this%mixing_ratios_conc=mixing_ratios_conc !> chapuza
    call this%allocate_mixing_waters_indices()
    do 
        read(unit,*) label
        if (label.eq.'end') then
            exit
        else if (label.eq.'MIXING RATIOS') then
            i=0 !> counter targets
            do
                i=i+1
                if (i<=this%spatial_discr%Num_targets) then
                    read(unit,*) tar_dim, (mixing_ratios_conc%cols(i)%col_1(j), j=1,tar_dim)
                    !print *, mixing_ratios_conc%cols(i)%col_1
                else
                    exit
                end if
            end do
            this%mixing_ratios_conc=mixing_ratios_conc !> chapuza
        else if (label.eq.'MIXING WATERS') then
            i=0 !> counter number of domain target waters
            do while (i<this%mixing_waters_indices%num_cols)
                i=i+1 !> update number of domain target waters
                !if (this%mixing_waters_indices%cols(i)%dim<this%mixing_waters_indices%num_cols-1) then
                !    read(unit,*) mix_wat_ind, (this%mixing_waters_indices%cols(i)%col_1(j), & !> mixing waters indices (we assume that the mixing waters are already ordered in the file)
                !        j=1,this%mixing_waters_indices%cols(i)%dim)
                !else
                    allocate(mix_wat_indices(this%mixing_ratios_conc%cols(i)%dim))
                    read(unit,*) mix_wat_ind, mix_wat_indices
                    this%mixing_waters_indices%cols(I)%col_1(1)=mix_wat_ind !> current target water
                    do j=1,this%mixing_ratios_conc%cols(i)%dim
                        if (mix_wat_indices(j).eq.mix_wat_ind) then
                            pos=j !> position of the current mixing water in the mixing ratios vector
                            this%mixing_ratios_conc%cols(I)%col_1(1)=mixing_ratios_conc%cols(I)%col_1(J) !> mixing ratio of the current target water
                            this%mixing_waters_indices%cols(I)%col_1(this%mixing_waters_indices%cols(I)%dim-1)=pos-1 !> number of upstream mixing waters
                            this%mixing_waters_indices%cols(I)%col_1(this%mixing_waters_indices%cols(I)%dim)=&
                                this%mixing_ratios_conc%cols(i)%dim-pos !> number of downstream mixing waters
                            exit
                        end if
                    end do
                    do j=1,this%mixing_ratios_conc%cols(i)%dim
                        if (j<pos) then !> upstream mixing waters
                            this%mixing_ratios_conc%cols(I)%col_1(1+j)=mixing_ratios_conc%cols(I)%col_1(J)
                            this%mixing_waters_indices%cols(I)%col_1(1+j)=mix_wat_indices(J)
                        else if (j>pos) then !> downstream mixing waters
                            this%mixing_ratios_conc%cols(I)%col_1(j)=mixing_ratios_conc%cols(I)%col_1(J)
                            this%mixing_waters_indices%cols(I)%col_1(j)=mix_wat_indices(J)
                        end if
                    end do
                    deallocate(mix_wat_indices)
                !end if
            end do
        else
            continue
        end if
    end do
    close(unit)
!> Chapuza
    !if (this%time_discr%int_method.eq.1) then !> Euler explicit
    !    do i=1,this%mixing_ratios_Rk%num_cols
    !        this%mixing_ratios_Rk%cols(i)%col_1=1d0
    !        print *, this%mixing_ratios_Rk%cols(i)%col_1
    !    end do
    !else
    !Delta_t=this%time_discr%get_Delta_t() !> we get time step (we assume uniform)
    if (this%mixing_ratios_conc%num_cols>=this%spatial_discr%Num_targets) then !> Euler implicit for transport
        do i=1,this%mixing_ratios_Rk%num_cols
            this%mixing_ratios_Rk%cols(i)%col_1=this%mixing_ratios_conc%cols(i)%col_1 !> we assume Euler fully implicit
            this%mixing_ratios_Rk_init%cols(i)%col_1=this%mixing_ratios_Rk%cols(i)%col_1
            !print *, this%mixing_ratios_Rk%cols(i)%dim
        end do
    else !> Euler explicit for transport
        do i=1,this%mixing_ratios_Rk%num_cols
            this%mixing_ratios_Rk%cols(i)%col_1=1d0
            this%mixing_ratios_Rk_init%cols(i)%col_1=this%mixing_ratios_Rk%cols(i)%col_1
            !print *, this%mixing_ratios_Rk%cols(i)%dim
        end do
    end if
end subroutine