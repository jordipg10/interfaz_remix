!> Writes data and final results of reactive transport problem
!> Transport data is written to {root}_transport.out
!> Chemistry data and results are written to {root}_chemistry.out
subroutine write_RT(this,dir,root,path_py)
    use RT_m, only: RT_c, RT_1D_transient_c, RT_2D_transient_c
    implicit none
    class(RT_c), intent(in) :: this                      !> reactive transport object
    character(len=*), intent(in) :: dir                      !> directory of output file
    character(len=*), intent(in) :: root                     !> root of output file
    character(len=*), intent(in), optional :: path_py        !> path for Python files

    integer(kind=4) :: unit,i

    select type (this)
    type is (RT_1D_transient_c)
        !> Write transport data to {root}_transport.out
        call this%transport%write_transport_data_WMA(dir, root)
        !> Re-open transport file in append mode for additional transport data
        open(newunit=unit, file=dir//root//'_transport.out', status='old', &
             form='formatted', position='append')
        if (this%transport%spatial_discr%num_targets<=200) then
            !> Write Lagrangian particle positions
            if (this%transport%Lagr_flag) then
                write(unit,"(2x,'Initial water particles positions:'/)")  !> write initial Lagrangian particle positions header
                do i=1,this%chemistry%num_target_waters_init !> loop over initial target waters
                    write(unit,"(2x,I5,ES15.5)") this%target_waters_init(i)%id, &
                        this%target_waters_init(i)%pos(1) !> write particle ID and spatial position [L]
                end do
                write(unit,"(/,2x,'Final water particles positions:'/)")  !> write final Lagrangian particle positions header
                do i=1,this%chemistry%num_target_waters !> loop over final target waters
                    write(unit,"(2x,I5,ES15.5)") this%target_waters(i)%id, &
                        this%target_waters(i)%pos(1) !> write particle ID and final position [L]
                end do
            end if
        end if
        close(unit)
        !> Write chemistry data and results to {root}_chemistry.out
        call this%chemistry%write_chemistry(dir, root)
    type is (RT_2D_transient_c)
        !> Write transport data to {root}_transport.out
        call this%transport%write_transport_data_WMA(dir, root)
        !> Re-open transport file in append mode for additional transport data
        open(newunit=unit, file=dir//root//'_transport.out', status='old', &
             form='formatted', position='append')
        if (this%transport%spatial_discr%num_targets<=200) then
            !> Write Lagrangian particle positions
            if (this%transport%Lagr_flag) then
                write(unit,"(2x,'Initial water particles positions:'/)")
                do i=1,this%chemistry%num_target_waters_init
                    write(unit,"(2x,I5,2ES15.5)") this%target_waters_init(i)%id, &
                        this%target_waters_init(i)%pos(1), this%target_waters_init(i)%pos(2)
                end do
                write(unit,"(/,2x,'Final water particles positions:'/)")
                do i=1,this%chemistry%num_target_waters
                    write(unit,"(2x,I5,2ES15.5)") this%target_waters(i)%id, &
                        this%target_waters(i)%pos(1), this%target_waters(i)%pos(2)
                end do
            end if
        end if
        close(unit)
        !> Write chemistry data and results to {root}_chemistry.out
        call this%chemistry%write_chemistry(dir, root)
    end select
end subroutine write_RT