!> Reads 'kinetics_modif.dat' database
subroutine read_kinetics_DB(this,path,unit)
    use chem_system_m, only: chem_system_c, mineral_c, kin_mineral_c, kin_mineral_params_c, aq_species_c
    implicit none
    class(chem_system_c) :: this                    !> chemical system
    character(len=*), intent(in) :: path            !> path to database
    integer(kind=4), intent(in) :: unit             !> file unit
    
    integer(kind=4) :: i,j,k,l,N_k,N_c,mineral_ind,cat_ind
    integer(kind=4), allocatable :: global_stoich_counter(:),global_stoich_indices(:),aux_ind(:)
    logical :: flag
    real(kind=8) :: unk,logK,valence,act_E,supersat_thr,cst
    real(kind=8), allocatable :: global_stoich_coeffs(:),aux_coeffs(:)
    character(len=256) :: name,comments,label,filename
    character(len=256), allocatable :: str1(:),str2(:),species_str(:),global_stoich_names(:)
    
    type(mineral_c) :: mineral
    type(kin_mineral_c), allocatable :: min_kin_reacts(:)
    type(kin_mineral_params_c) :: min_params
    
    type(aq_species_c), allocatable :: catalysers(:)
    
    allocate(min_kin_reacts(this%num_minerals_kin))
    
    filename=trim(path)//'\kinetics_modif.dat'

    open(unit,file=filename,status='old',action='read')
    i=1 !> counter kinetic minerals chemical system
    j=1 !> counter number parallel reactions
    k=1 !> counter number catalysts
    read(unit,*) label
    do
        read(unit,*) mineral%name, N_k, N_c, act_E, supersat_thr, comments
        if (name=='null') then
            exit
        else
            call this%is_mineral_in_chem_syst(mineral,flag,mineral_ind)
            if (flag.eqv..true. .and. mineral_ind<=this%num_minerals_kin) then
                min_params%num_par_reacts=N_k
                min_params%num_cat=N_c
                min_params%act_energy=act_E*1484 !> reads kcal/mol so we have to convert to Joules
                min_params%supersat_threshold=supersat_thr
                call min_params%allocate_constants()
                call min_params%allocate_cat_indices()
                allocate(catalysers(min_params%num_cat))
                do j=1,min_params%num_par_reacts
                    !read(unit,*) min_params%k(j), min_params%theta(j), min_params%eta(j), min_params%num_cat, ((min_params%catalysers(l)%name, min_params%p(j,l)), l=1,min_params%num_cat)
                    read(unit,*) min_params%k(j), min_params%theta(j), min_params%eta(j)
                    do k=1,min_params%num_cat
                        read(unit,*) catalysers(k)%name, min_params%p(j,k)
                    end do
                    do k=1,min_params%num_cat
                        call this%aq_phase%is_species_in_aq_phase(catalysers(k),flag,cat_ind)
                        if (flag.eqv..false.) then
                            error stop 'Catalyser not found in chemical system'
                        else
                            min_params%cat_indices(k)=cat_ind
                        end if
                    end do
                end do
                call this%min_kin_reacts(mineral_ind)%set_mineral_params(min_params)
                if (i<this%num_minerals_kin) then
                    i=i+1
                else
                    exit
                end if
                deallocate(min_params%k,min_params%p,min_params%theta,min_params%eta,min_params%cat_indices,catalysers)
            else
                read(unit,*) cst
            end if
        end if
    end do
    close(unit)
end subroutine