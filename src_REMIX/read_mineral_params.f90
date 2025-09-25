!> Reads parameters for mineral kinetic reaction
subroutine read_mineral_params(this,react_name,filename)
    use kin_mineral_params_m
    implicit none
    class(kin_mineral_params_c) :: this
    character(len=*), intent(in) :: react_name
    character(len=*), intent(in) :: filename
            
    integer(kind=4) :: i,j,n_t,index,water_flag,n_paths
    real(kind=8) :: k
    real(kind=8), allocatable :: k_vec(:),k_inh(:),k_M(:)
    character(len=256) :: str,str2,str3,cat
    integer(kind=4), allocatable :: global_stoich_counter(:),global_stoich_indices(:),aux_ind(:)
    character(len=256), allocatable :: str_array(:),cat_array(:),global_stoich_names(:)
    real(kind=8), allocatable :: global_stoich_coeffs(:),aux_coeffs(:)
    logical :: flag
    type(species_c), allocatable :: inhibitors(:),TEAs(:) 
    
    error stop "Subroutine 'read_mineral_params' not implemented yet"
end subroutine