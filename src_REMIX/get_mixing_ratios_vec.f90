!> Gets the mixing ratios associated to a target
subroutine get_mixing_ratios_vec(mixing_ratios,j,mixing_ratios_vec)
    use chemistry_m
    implicit none
    class(real_array_c), intent(in) :: mixing_ratios
    !real(kind=8), intent(in) :: mixing_ratios(:) ! dim = n x 4 ([0,sub],diag,[super,0],revalence) if transport explicit or n x n if transport implicit
    integer(kind=4), intent(in) :: j ! target index
    real(kind=8), intent(out) :: mixing_ratios_vec(:)   ! dim=3 or 4
    
    mixing_ratios_vec=mixing_ratios%cols(j)%col_1
        !if (j==1) then
        !    !mixing_ratios_vec=[mixing_ratios%diag(j),mixing_ratios%super(j)]
        !    mixing_ratios_vec=[mixing_ratios(j,2),mixing_ratios(j,3),mixing_ratios(j,4)]
        !else if (j==size(mixing_ratios,1)) then
        !    !mixing_ratios_vec=[mixing_ratios%diag(j),mixing_ratios%sub(j-1)]
        !    mixing_ratios_vec=[mixing_ratios(j,2),mixing_ratios(j,1),mixing_ratios(j,4)]
        !else
        !    !mixing_ratios_vec=[mixing_ratios%diag(j),mixing_ratios%sub(j-1),mixing_ratios%super(j)]
        !    mixing_ratios_vec=[mixing_ratios(j,2),mixing_ratios(j,1),mixing_ratios(j,3),mixing_ratios(j,4)]
        !end if
end subroutine